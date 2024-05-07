import pandas as pd
import numpy as np
import os
import os.path
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError
from abc_lite.tools import *
import linecache
import traceback
import time
import pyranges as pr

### minimize the necessary codebase

pd.options.display.max_colwidth = 10000 #seems to be necessary for pandas to read long file names... strange

def get_tss_for_bed(bed):
    assert_bed3(bed)
    tss = bed['start'].copy()
    tss.loc[bed.loc[:,'strand'] == "-"] = bed.loc[bed.loc[:,'strand'] == "-",'end']

    return tss

def assert_bed3(df):
    assert(type(df).__name__ == "DataFrame")
    assert('chr' in df.columns)
    assert('start' in df.columns)
    assert('end' in df.columns)
    assert('strand' in df.columns)

#Kristy's version
def assign_enhancer_classes(enhancers, genes, tss_slop=500):

    # build pyranges df 
    tss_pyranges = df_to_pyranges(genes, start_col='tss', end_col='tss', start_slop=tss_slop, end_slop=tss_slop)
    gene_pyranges = df_to_pyranges(genes)

    def get_class_pyranges(enhancers, tss_pyranges = tss_pyranges, gene_pyranges = gene_pyranges): 
        '''
        Takes in PyRanges objects : Enhancers, tss_pyranges, gene_pyranges
        Returns dataframe with  uid (representing enhancer) and symbol of the gene/promoter that is overlapped'''

        #genes
        genic_enh = enhancers.join(gene_pyranges, suffix="_genic")
        genic_enh = genic_enh.df[['symbol','uid']].groupby('uid',as_index=False).aggregate(lambda x: ','.join(list(set(x))))
        
        #promoters
        promoter_enh = enhancers.join(tss_pyranges, suffix="_promoter")
        promoter_enh = promoter_enh.df[['symbol','uid']].groupby('uid',as_index=False).aggregate(lambda x: ','.join(list(set(x))))
        
        return genic_enh, promoter_enh

    # label everything as intergenic
    enhancers["class"] = "intergenic"
    enhancers['uid'] = range(enhancers.shape[0])
    enh = df_to_pyranges(enhancers)
 
    genes, promoters = get_class_pyranges(enh)
    enhancers = enh.df.drop(['Chromosome','Start','End'], axis=1)
    enhancers.loc[enhancers['uid'].isin(genes.uid), 'class'] = 'genic'
    enhancers.loc[enhancers['uid'].isin(promoters.uid), 'class'] = 'promoter' 
    
    enhancers["isPromoterElement"] = enhancers["class"] == "promoter"
    enhancers["isGenicElement"] = enhancers["class"] == "genic"
    enhancers["isIntergenicElement"] = enhancers["class"] == "intergenic"
  
    # Output stats
    print("Total enhancers: {}".format(len(enhancers)))
    print("         Promoters: {}".format(sum(enhancers['isPromoterElement'])))
    print("         Genic: {}".format(sum(enhancers['isGenicElement'])))
    print("         Intergenic: {}".format(sum(enhancers['isIntergenicElement'])))

    #Add promoter/genic symbol
    enhancers = enhancers.merge(promoters.rename(columns={'symbol':'promoterSymbol'}), on='uid', how = 'left').fillna(value={'promoterSymbol':""})
    enhancers = enhancers.merge(genes.rename(columns={'symbol':'genicSymbol'}), on='uid', how = 'left').fillna(value={'genicSymbol':""})
    enhancers.drop(['uid'], axis=1, inplace=True)

    # just to keep things consistent with original code 
    enhancers["name"] = enhancers.apply(lambda e: "{}|{}:{}-{}".format(e["class"], e.chr, e.start, e.end), axis=1)
    return enhancers


#
bed_extra_colnames = ["name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]
#JN: 9/13/19: Don't assume chromosomes start with 'chr'
#chromosomes = ['chr' + str(entry) for entry in list(range(1,23)) + ['M','X','Y']]   # should pass this in as an input file to specify chromosome order
def read_bed(filename, extra_colnames=bed_extra_colnames, chr=None, sort=False, skip_chr_sorting=True):
    skip = 1 if ("track" in open(filename, "r").readline()) else 0
    names = ["chr", "start", "end"] + extra_colnames
    result = pd.read_table(filename, names=names, header=None, skiprows=skip, comment='#')
    result = result.dropna(axis=1, how='all')  # drop empty columns
    assert result.columns[0] == "chr"

    #result['chr'] = pd.Categorical(result['chr'], chromosomes, ordered=True)
    result['chr'] = pd.Categorical(result['chr'], ordered=True)
    if chr is not None:
        result = result[result.chr == chr]
    if not skip_chr_sorting:
        result.sort_values("chr", inplace=True)
    if sort:
        result.sort_values(["chr", "start", "end"], inplace=True)
    return result


def my_qnorm(df, separate_promoters=True):
    df['qnormed'] = df['activity_base']
    if not separate_promoters:
        df['qnormed'] = df['activity_base'].rank()/df.shape[0]
    else:
        this_idx = df.index[np.logical_or(df['class'] == "tss", df['class'] == "promoter")]
        df.loc[this_idx, 'qnormed'] = df.loc[this_idx, 'activity_base'].rank()/len(this_idx)
        this_idx = df.index[np.logical_or(df['class'] == "tss", df['class'] == "promoter")]
        df.loc[this_idx, 'qnormed'] = df.loc[this_idx, 'activity_base'].rank()/len(this_idx)
    df['un_normed'] = df['activity_base']
    df['activity_base'] = df['qnormed']
    return df

def run_qnorm(df, qnorm, qnorm_method = "rank", separate_promoters = True):
    # Quantile normalize epigenetic data to a reference
    #
    # Option to qnorm promoters and nonpromoters separately

    if qnorm is None:
        if 'H3K27ac.RPM' in df.columns: df['normalized_h3K27ac'] = df['H3K27ac.RPM']
        if 'DHS.RPM' in df.columns: df['normalized_dhs'] = df['DHS.RPM']
        if 'ATAC.RPM' in df.columns: df['normalized_atac'] = df['ATAC.RPM']
    else:
        qnorm = pd.read_csv(qnorm, sep = "\t")
        nRegions = df.shape[0] 
        col_dict = {'DHS.RPM' : 'normalized_dhs', 'ATAC.RPM' : 'normalized_atac', 'H3K27ac.RPM' : 'normalized_h3K27ac'}

        for col in set(df.columns & col_dict.keys()):
            #if there is no ATAC.RPM in the qnorm file, but there is ATAC.RPM in enhancers, then qnorm ATAC to DHS
            if col == 'ATAC.RPM' and 'ATAC.RPM' not in qnorm.columns:
                qnorm['ATAC.RPM'] = qnorm['DHS.RPM']

            if not separate_promoters:
                qnorm = qnorm.loc[qnorm['enh_class' == "any"]]
                if qnorm_method == "rank":
                    interpfunc = interpolate.interp1d(qnorm['rank'], qnorm[col], kind='linear', fill_value='extrapolate')
                    df[col_dict[col]] = interpfunc((1 - df[col + ".quantile"]) * nRegions).clip(0)
                elif qnorm_method == "quantile":
                    interpfunc = interpolate.interp1d(qnorm['quantile'], qnorm[col], kind='linear', fill_value='extrapolate')
                    df[col_dict[col]] = interpfunc(df[col + ".quantile"]).clip(0)
            else:
                for enh_class in ['promoter','nonpromoter']:
                    this_qnorm = qnorm.loc[qnorm['enh_class'] == enh_class]

                    #Need to recompute quantiles within each class
                    if enh_class == 'promoter':
                        this_idx = df.index[np.logical_or(df['class'] == "tss", df['class'] == "promoter")]
                    else:
                        this_idx = df.index[np.logical_and(df['class'] != "tss" , df['class'] != "promoter")]
                    df.loc[this_idx, col + enh_class + ".quantile"] = df.loc[this_idx, col].rank()/len(this_idx)

                    if qnorm_method == "rank":
                        interpfunc = interpolate.interp1d(this_qnorm['rank'], this_qnorm[col], kind='linear', fill_value='extrapolate')
                        df.loc[this_idx, col_dict[col]] = interpfunc((1 - df.loc[this_idx, col + enh_class + ".quantile"]) * len(this_idx)).clip(0)
                    elif qnorm_method == "quantile":
                        interpfunc = interpolate.interp1d(this_qnorm['quantile'], this_qnorm[col], kind='linear', fill_value='extrapolate')
                        df.loc[this_idx, col_dict[col]] = interpfunc(df.loc[this_idx, col + enh_class + ".quantile"]).clip(0)

    return df
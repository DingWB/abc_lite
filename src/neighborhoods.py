import pandas as pd
import numpy as np
from scipy import interpolate
import pysam
import os
import os.path
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError
from tools import *
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

    # import pdb
    # pdb.Pdb(stdout=sys.__stdout__).set_trace()
    # pdb.set_trace()

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


def count_single_feature_for_bed(df, bed_file, genome_sizes, feature_bam, feature, directory, filebase, skip_rpkm_quantile, force, use_fast_count):
    orig_shape = df.shape[0]
    feature_name = feature + "." + os.path.basename(feature_bam)
    feature_outfile = os.path.join(directory, "{}.{}.CountReads.bedgraph".format(filebase, feature_name))

    if force or (not os.path.exists(feature_outfile)) or (os.path.getsize(feature_outfile) == 0):
        print("Regenerating", feature_outfile)
        print("Counting coverage for {}".format(filebase + "." + feature_name))
        run_count_reads(feature_bam, feature_outfile, bed_file, genome_sizes, use_fast_count)
    else:
        print("Loading coverage from pre-calculated file for {}".format(filebase + "." + feature_name))

    domain_counts = read_bed(feature_outfile)
    score_column = domain_counts.columns[-1]

    total_counts = count_total(feature_bam)

    domain_counts = domain_counts[['chr', 'start', 'end', score_column]]
    featurecount = feature_name + ".readCount"
    domain_counts.rename(columns={score_column: featurecount}, inplace=True)
    domain_counts['chr'] = domain_counts['chr'].astype('str')

    df = df.merge(domain_counts.drop_duplicates())
    #df = smart_merge(df, domain_counts.drop_duplicates())

    assert df.shape[0] == orig_shape, "Dimension mismatch"

    df[feature_name + ".RPM"] = 1e6 * df[featurecount] / float(total_counts)

    if not skip_rpkm_quantile:
        df[featurecount + ".quantile"] = df[featurecount].rank() / float(len(df))
        df[feature_name + ".RPM.quantile"] = df[feature_name + ".RPM"].rank() / float(len(df))
        df[feature_name + ".RPKM"] = 1e3 * df[feature_name + ".RPM"] / (df.end - df.start).astype(float)
        df[feature_name + ".RPKM.quantile"] = df[feature_name + ".RPKM"].rank() / float(len(df))

    return df[~ df.duplicated()]

def average_features(df, feature, feature_bam_list, skip_rpkm_quantile):
    feature_RPM_cols = [feature + "." + os.path.basename(feature_bam) + '.RPM' for feature_bam in feature_bam_list]

    df[feature + '.RPM'] = df[feature_RPM_cols].mean(axis = 1)
    
    if not skip_rpkm_quantile:
        feature_RPKM_cols = [feature + "." + os.path.basename(feature_bam) + '.RPKM' for feature_bam in feature_bam_list]
        df[feature + '.RPM.quantile'] = df[feature + '.RPM'].rank() / float(len(df))
        df[feature + '.RPKM'] = df[feature_RPKM_cols].mean(axis = 1)
        df[feature + '.RPKM.quantile'] = df[feature + '.RPKM'].rank() / float(len(df))

    return df

# From /seq/lincRNA/Jesse/bin/scripts/JuicerUtilities.R
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


def read_bedgraph(filename):
    read_bed(filename, extra_colnames=["score"], skip_chr_sorting=True)

def count_bam_mapped(bam_file):
    # Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
    # chromosomes = ['chr' + str(x) for x in range(1,23)] + ['chrX'] + ['chrY']
    command = ("samtools idxstats " + bam_file)
    data = check_output(command, shell=True)
    lines = data.decode("ascii").split("\n")
    #vals = list(int(l.split("\t")[2]) for l in lines[:-1] if l.split("\t")[0] in chromosomes)
    vals = list(int(l.split("\t")[2]) for l in lines[:-1])
    if not sum(vals) > 0:
        raise ValueError("Error counting BAM file: count <= 0")
    return sum(vals)

def count_tagalign_total(tagalign):
    #result = int(check_output("zcat " + tagalign + " | wc -l", shell=True))
    result = int(check_output("zcat {} | grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' | wc -l".format(tagalign), shell=True))
    assert (result > 0)
    return result

def count_bigwig_total(bw_file):
    from pyBigWig import open as open_bigwig
    bw = open_bigwig(bw_file)
    result = sum(l * bw.stats(ch, 0, l, "mean")[0] for ch, l in bw.chroms().items())
    assert (abs(result) > 0)  ## BigWig could have negative values, e.g. the negative-strand GroCAP bigwigs
    return result

def count_total(infile):
    if infile.endswith(".tagAlign.gz") or infile.endswith(".tagAlign.bgz"):
        total_counts = count_tagalign_total(infile)
    elif infile.endswith(".bam"):
        total_counts = count_bam_mapped(infile)
    elif isBigWigFile(infile):
        total_counts = count_bigwig_total(infile)
    else:
        raise RuntimeError("Did not recognize file format of: " + infile)

    return total_counts

def parse_params_file(args):
    # Parse parameters file and return params dictionary
    params = {}

    params["default_accessibility_feature"] = determine_accessibility_feature(args)
    params["features"] = get_features(args)

    if args.expression_table:
        params["expression_table"] = args.expression_table.split(",")
    else:
        params["expression_table"] = ''

    return(params)

def get_features(args):
    features = {}

    if args.H3K27ac:
        features['H3K27ac'] = args.H3K27ac.split(",")
    
    if args.ATAC:
        features['ATAC'] = args.ATAC.split(",")
    
    if args.DHS:
        features['DHS'] = args.DHS.split(",")

    if args.supplementary_features is not None:
        supp = pd.read_csv(args.supplementary_features, sep="\t")
        for idx,row in supp.iterrows():
            features[row['feature_name']] = row['file'].split(",")

    return features

def determine_accessibility_feature(args):
    if args.default_accessibility_feature is not None:
        return args.default_accessibility_feature
    elif (not args.ATAC) and (not args.DHS):
        raise RuntimeError("Both DHS and ATAC have been provided. Must set one file to be the default accessibility feature!")
    elif args.ATAC:
        return "ATAC"
    elif args.DHS:
        return "DHS"
    else:
        raise RuntimeError("At least one of ATAC or DHS must be provided!")

def compute_activity(df, access_col):
    if access_col == "DHS":
        if 'H3K27ac.RPM' in df.columns:
            df['activity_base'] = np.sqrt(df['normalized_h3K27ac'] * df['normalized_dhs'])
            df['activity_base_no_qnorm'] = np.sqrt(df['H3K27ac.RPM'] * df['DHS.RPM'])
        else:
            df['activity_base'] = df['normalized_dhs']
            df['activity_base_no_qnorm'] = df['DHS.RPM']
    elif access_col == "ATAC":
        if 'H3K27ac.RPM' in df.columns:
            df['activity_base'] = np.sqrt(df['normalized_h3K27ac'] * df['normalized_atac'])
            df['activity_base_no_qnorm'] = np.sqrt(df['H3K27ac.RPM'] * df['ATAC.RPM'])
        else:
            df['activity_base'] = df['normalized_atac']
            df['activity_base_no_qnorm'] = df['ATAC.RPM']
    else:
        raise RuntimeError("At least one of ATAC or DHS must be provided!")

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

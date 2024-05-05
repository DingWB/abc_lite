#!python
import os
import sys
# sys.path.insert(0, '/data/earmand/projects/dp_hic_revision/abc_model/ABC-Enhancer-Gene-Prediction/src')
from abc_lite.neighborhoods import assign_enhancer_classes, get_tss_for_bed
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def sanitize_cols(df):
    df = df.loc[:,~df.columns.duplicated()].copy()
    return df

def format_gene_table(in_path, outdir, chrom_sizes,
                     header=False, tss_slop=1000):
    '''
    Read in gene table and format it for ABC model
    
    inputs:
        in_path - the input path of the gene table
        header  - whether the input file has a header
    '''
    # needed columns
    # isExpressed is determined during prediction. assuming we have expression and
    # don't calculate activity quantiles, we will put activity qunatiles to zero
    # the means we just need the chr, symbol, tss, and expression columns
    # ['chr','symbol','tss','Expression','PromoterActivityQuantile','isExpressed']

    if header:
        gene_table = pd.read_csv(in_path, sep='\t')
        if 'chrom' in gene_table.columns:
            gene_table.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    else:
        gene_table = pd.read_csv(in_path, sep='\t', header=None)
        gene_table.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    
    tss = get_tss_for_bed(gene_table)

    gene_table['tss'] = tss
    save_genes = gene_table.loc[:,['chr',  'start', 'end', 'name', 'start', 'score', 'strand',  'tss']]
    save_genes['symbol'] = save_genes['name']
    save_genes['Expression'] = save_genes['score']
    save_genes['PromoterActivityQuantile'] = 0
    save_genes = sanitize_cols(save_genes)
    save_genes.to_csv(os.path.join(outdir, "GeneList.txt"), sep='\t', index=False, header=True, float_format="%.6f")
    # save_genes = gene_table.loc[:,['chr', 'name', 'start', 'end', 'score']]
    return save_genes

# annotate_genes_with_features
# genes['Expression.quantile'] = genes['Expression'].rank(method='average', na_option="top", ascending=True, pct=True)
def format_enhancer_table(in_path, gene_table, outdir, 
                        qnorm=False, header=False, 
                        signal_col=4, tss_slop=1000) -> None:
    # needed columns
    # ['chr','start','end','name','class','activity_base']
    '''
    Format enhancer table to be used in ABC model
    
    inputs:
        in_path      - the input path of the bed format enhancer activity quantification
        gene_table   - dataframe of genes
        qnorm        - whether to quantile normalize the activity values
        header       - whether the input file has a header
        signal_col - the name of the column with the ehancer singal
    '''
    # read in the enhancer table
    if header==False:
        enhancer_table = pd.read_csv(in_path, sep='\t', header=None)
        enhancer_table = enhancer_table.iloc[:,[0,1,2,signal_col]]
        enhancer_table.columns = ['chr', 'start', 'end', signal_col]
    else:
        enhancer_table = pd.read_csv(in_path, sep='\t')
        # try:
        if 'chrom' in enhancer_table.columns:
            enhancer_table = enhancer_table.loc[:,['chrom', 'start', 'end', signal_col]]
            enhancer_table.columns = ['chr', 'start', 'end', signal_col]
    enhancer_table['activity_base'] = enhancer_table[signal_col].values

    # quantile normalize the activity values

    #     enhancer_table['activity_base'] = enhancer_table['activity_base'].rank() / float(len(enhacer_table))
    # santatize the enhancer table before pyranges
    enhancer_table = sanitize_cols(enhancer_table)
    enhancer_table = assign_enhancer_classes(enhancer_table,
                                         gene_table, 
                                         tss_slop)
    if qnorm:
    # todo implement quantile normalization
        pass
    enhancer_table.to_csv(os.path.join(outdir, "EnhancerList.txt"),
                sep='\t', index=False, header=True, float_format="%.6f")


def main():
    parser = ArgumentParser()
    parser.add_argument('--atac', type=str,
         required=True, help='Path to ATAC-seq data in .bed format')
    parser.add_argument('--rna', type=str,
            required=True, help='Path to RNA-seq data in 6 column bed format')
    parser.add_argument('--hic', type=str,
            required=True, help='Path to Hi-C data in .hic format')
    parser.add_argument('--out', type=str,
            required=True, help='Output directory')
    parser.add_argument('--chrom_sizes', type=str,
            required=True, help='Path to chrom sizes file')
    parser.add_argument('--hic_resolution', type=int,
            default=10000, help='Resolution of Hi-C data')
    parser.add_argument('--tss_slop', type=int,
            default=1000, help='Distance around tss to assign peaks as tss')
    parser.add_argument('--header', action='store_true',
            help='Whether the input files have headers',
             default=False)
    parser.add_argument('--qnorm', action='store_true',
            help='Wether to qunatile normalize atac signal')
    parser.add_argument('--signal_col', 
            default=4, help='Column with signal values in atac file')
    parser.add_argument('--juicebox_path', type=str, default='$JUICERTOOLS',
            help='Path to juicebox jar file')
    parser.add_argument('--src_path', type=str,
            default='/data/earmand/projects/dp_hic_revision/abc_model/ABC-Enhancer-Gene-Prediction/src')
    parser.add_argument('--window_size', type=int, default=5000000,
            help='Window size around gene to look for enhancers')
    parser.add_argument('--signal_column', default=4,
            help='Column with signal values in atac file')
    parser.add_argument('--gene_quantile', type=float, default=0.4,
            help='the minimum expression qunatile for which abc links are calculated')
    parser.add_argument('--threshold', type=float, default=0.2,
            help='Threshold for ABC model peak-gene links')
    args = parser.parse_args()

    # read in gene table

    gene_table = format_gene_table(args.rna, args.out, args.chrom_sizes,
                                    args.header, args.tss_slop)

    # read in enhancer table
    format_enhancer_table(in_path = args.atac, gene_table = gene_table, 
                qnorm=args.qnorm,
             header=args.header, signal_col=args.signal_column, 
             tss_slop =args.tss_slop, outdir=args.out)

    # dump hic
    os.makedirs(f'{args.out}/hic_dump', exist_ok=True)
    hic_command = ( f'python {args.src_path}/juicebox_dump.py '
     f'--hic_file {args.hic} --resolution {args.hic_resolution} '
     f'--outdir {args.out}/hic_dump --juicebox "java -jar {args.juicebox_path}"'
    )
    print('hic command: \n', hic_command)
    os.system(hic_command)

    # run abc
    abc_command = (f'python {args.src_path}/predict.py ' 
                   f'--genes {args.out}/GeneList.txt '
                   f'--enhancers {args.out}/EnhancerList.txt '
                   f'--threshold 0.2 '
                   f'--HiCdir {args.out}/hic_dump --outdir {args.out} '
                   f'--window {args.window_size} '
                   '--make_all_putative '
                   f'--promoter_activity_quantile_cutoff {args.gene_quantile} '
                   f'--hic_resolution {args.hic_resolution} '

    )
    print('abc command: \n', abc_command)
    os.system(abc_command)


if __name__ == '__main__':
    main()







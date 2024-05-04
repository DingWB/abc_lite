# ABC (activity by contact model) lite

This repo implements a lightweight version of the abc model built on the assumption of users having a peak set/feature table. 

## limitations
Currently this implementation does not support quantile normalization or multifeature activity (e.g.  both H3k27ac & ATAC). 
It operates under the assumption you have a gene expression measure, epigenome measure, and HiC

#### unimplemented features
- Estimate HiC from genomic distance
- quantile normalization

## Installation
Installation relies on pip and conda.

1. clone the repo
2. setup conda env
suggested python 3.9, not tested with other versions
```
conda create -n abc_lite python=3.9
conda activate abc_lite
```
3. install requirements
this additionally downloads juicertools and sets the environmental variable JUICERTOOLS in the environment
```
bash ./conda_install.sh
```
4. Finish setup with pip
```
python3 -m pip install .
```

## citation
If you use the ABC model in published research, please cite:

[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

[2] Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F, Mualim K, Natri HM, Weeks EM, Munson G, Kane M, Kang HY, Cui A, Ray JP, Eisenhaure TM, Collins RL, Dey K, Pfister H, Price AL, Epstein CB, Kundaje A, Xavier RJ, Daly MJ, Huang H, Finucane HK, Hacohen N, Lander ES, Engreitz JM. Genome-wide enhancer maps link risk variants to disease genes. Nature. 2021 May;593(7858):238-243. doi: 10.1038/s41586-021-03446-x

## Requirements
For each cell-type, the inputs to the ABC model are:


 * Required Inputs
 	* Activity bedfile (e.g. ATAC, H3k27ac)  
		* required columns "chr", "start", "end" and a signal column (asummed rpkm/cpm)
		* header optional, if included should match requried columns
	* Gene bedfile bed6 format bed file
		* requried columns "chr", "start", "end", "name" "signal" "strand" (signal assumed normalized expression )
		* header optional, if included should match reqired columns
 	* Hi-C data 
		*Hi-C data in juicer (.hic) format

### Dependencies

Working on reformatting dependancies.
'''
Python (3.9)
'''


### Example usage

```
usage: run_abc_light.py [-h] --atac ATAC --rna RNA --hic HIC --out OUT --chrom_sizes CHROM_SIZES [--hic_resolution HIC_RESOLUTION] [--tss_slop TSS_SLOP] [--header] [--qnorm] [--signal_col SIGNAL_COL] [--juicebox_path JUICEBOX_PATH] [--src_path SRC_PATH] [--window_size WINDOW_SIZE] [--signal_column SIGNAL_COLUMN]
                        [--gene_quantile GENE_QUANTILE] [--threshold THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  --atac ATAC           Path to ATAC-seq data in .bed format
  --rna RNA             Path to RNA-seq data in 6 column bed format
  --hic HIC             Path to Hi-C data in .hic format
  --out OUT             Output directory
  --chrom_sizes CHROM_SIZES
                        Path to chrom sizes file
  --hic_resolution HIC_RESOLUTION
                        Resolution of Hi-C data
  --tss_slop TSS_SLOP   Distance around tss to assign peaks as tss
  --header              Whether the input files have headers
  --qnorm               Wether to qunatile normalize atac signal
  --signal_col SIGNAL_COL
                        Column with signal values in atac file
  --juicebox_path JUICEBOX_PATH
  --src_path SRC_PATH
  --window_size WINDOW_SIZE
                        Window size around gene to look for enhancers
  --signal_column SIGNAL_COLUMN
                        Column with signal values in atac file
  --gene_quantile GENE_QUANTILE
                        the minimum expression qunatile for which abc links are calculated
  --threshold THRESHOLD
                        Threshold for ABC model peak-gene links

```

## Description of the ABC Model

The Activity by Contact (ABC) model is designed to represent a mechanistic model in which enhancers activate gene transcription upon enhancer-promoter contact. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency with which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the contribution of a specific enhancer to a gene’s expression should depend on the surrounding context (ie, the strength and contact frequency of other enhancers for the gene).

To convert this conceptual framework into a practical score (which can be applied genome-wide), we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G /  Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Operationally, Activity (A) is defined as the geometric mean of the read counts of DNase-seq and H3K27ac ChIP-seq at an element E, and Contact (C) as the KR normalized Hi-C contact frequency between E and the promoter of gene G. Elements are defined as ~500bp regions centered on DHS peaks. 

Note that the ABC model only considers candidate elements and genes on the same chromosome. It does not make interchromosomal predictions.


* Accurate transcription start site annotations are critical. The ABC model uses the TSS of the gene in order to assign enhancer-promoter contact frequency. If the TSS annotation is inaccurate (off by >5kb) it will lead to inaccurate predictions.
* The size of candidate enhancer elements is important. For example, if two candidate regions are merged, then the ABC score of the merged region will be approximately the sum of the ABC scores for each individual region.
* In our testing the ABC model typically predicts on average ~3 distal enhancers per expressed gene. If you run the model on a cell type and find a large deviation from this number (say <2 or >5) this may mean the ABC model is not well calibrated in the cell type. Typical remedies are to use quantile normalization, scale Hi-C or to lower/raise the cutoff on the ABC score.

## Contact
Please submit a github issue with any questions or if you experience any issues/bugs. 



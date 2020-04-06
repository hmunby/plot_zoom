gemma2zoom 
====

gemma2zoom is a tool for visualising results from genome-wide analyses that test individual SNPs. The default parameters assume the results file has the same format as GEMMA output.  
It produces a manhattan plot of a local region around a focal SNP, colouring SNPs by their level of LD with the focal SNP. 

The script itself is a wrapper around region_plotter.R which does the plotting, this R script must remain in the same directory as gemma2zoom for it to work. 

Dependencies:

[PLINK v1.9](https://www.cog-genomics.org/plink2)

R packages:
Bioconductor
Rtracklayer


## Before you start 
____

Make sure that PLINK is installed and that the executable (plink) is in your path or the current working directory.  

Install Bioconductor & rtracklayer: 

``if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")`

`BiocManager::install("rtracklayer")`

### Generate PLINK BED files 

You will first need to produce PLINK BED files from the VCF which you would like LD values to be calculated from. 

This will generally be the same for each analysis that uses the same dataset; for my Lake Masoko analyses (all_CalMas) I have included these in the examples subdirectory. 

If using a different dataset then first use PLINK to generate these files: 

E.g. The example Lake Masoko files were generated as follows:
`plink --vcf dataset_name.vcf.gz --make-bed --out dataset_name`

`plink --vcf ~/rds/rds-rd109-durbin-group/projects/cichlid/massoko/all_CalMas/masked_snps/sex/all_CalMas_variants_genome_clean.vcf --make-bed --set-missing-var-ids @:# --out all_CalMas`

Or generally (input vcf can be gzipped or uncompressed):
`plink --vcf dataset_name.vcf.gz --make-bed --set-missing-var-ids @:# --out dataset_name`


## Usage
___

You must give the full path to the gemma2zoom executable script for now (i.e. don't put it in your path, it won't work)

### Options

Executing gemma2zoom without any options will print a list of options & descriptors:
Options required are denoted with `*` 
`* -i | --input      <string>        Path to input file containing results to be plotted.
  -h | --header     <0|1>           First line of input file is not (0) or is (1) a header. (default assumes header present)
* -t | --testcol    <int>           Number of column (1-indexed) containing test statistic to be plotted.
-x | --chrcol     <int>           Number of column (1-indexed) indicating site chromosome. (default is col 1)
-y | --poscol     <int>           Number of column (1-indexed) indicating site position. (default is col 3)
* -b | --plinkbed   <string>        Path to PLINK BED file (generated directly from VCF) to be used for LD calculations (prefix only).
* -a | --annot      <string>        Path to GTF annotation file.
-c | --fschr      <int>           Focal SNP: chromosome. (default is FIRST POSITION in assoc input file i.e. top SNP only if sorted by pval)
-p | --fspos      <int>           Focal SNP: position. (default is FIRST POSITION in assoc input file i.e. top SNP only if sorted by pval)
-w | --window     <int>           Size of window either side of focal SNP to be plotted in base pairs. (default is 50,000 bp)
-s | --sigline    <float|string>  Significance threshold value OR 'bfc' to plot Bonferroni corrected p=0.05. If unused no sig line plotted.
-l | --logtrans   <0|1>           Plot absolute statistic values (0) or their -log10 transform (1). (default is to log transform)
* -o | --outpref    <string>        Path/filename for output files (Plot and PLINK LD calculation files).`

### Examples

`./gemma2zoom.sh -i example/masoko_sex_GWAS.assoc.txt -t 14 -b example/all_CalMas -a example/Astatotilapia_calliptera.fAstCal1.2.99.chr.gtf -c 19 -p 21581905 -w 100000 -s bfc -o example/masoko_sex_GWAS_chr19`

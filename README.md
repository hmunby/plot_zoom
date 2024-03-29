plot_zoom
====

![Example output plot](example.png)

plot_zoom  is a tool for visualising results from any genome-wide analysis that produces values for individual sites. Default parameters assumes an input file in the format of a [GEMMA](https://github.com/genetics-statistics/GEMMA) association output file but parameters can be adjusted to accomodate variations to this format. 

plot_zoom outputs:
- Figure (png) of the region with 3 panels: 
	- A manhattan plot of a local region around a focal SNP, colouring SNPs by their level of LD with the focal SNP.
	- A plot showing the distribution of plotted SNPs.  
	- Gene annotations in the displayed region.
- List of coordinates of SNPs in the window and their r^2 LD values with the focal SNP 

NB: As the script itself is a wrapper around region_plotter.R which performs all plotting, this R script must remain in the same directory as the plot_zoom executable.

Dependencies:

[PLINK v1.9](https://www.cog-genomics.org/plink2)

R packages: 
[Bioconductor: rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)

Installation: 
`git pull https://github.com/hmunby/plot_zoom/`
Ready to use once repository has been pulled. Keep both .R and .sh scripts in the same directory and give path to the plot_zoom.sh to run. 

## Before you start
____

Ensure that PLINK is installed and that the executable (plink) is in your path or the current working directory.

Install Bioconductor & rtracklayer:

`if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")`

`BiocManager::install("rtracklayer")`

### Generate PLINK BED files

You will first need to produce PLINK BED files from the VCF which you would like LD values to be calculated from.

This will generally be the same for each analysis performed on the same dataset; for my Lake Masoko Astatotilapia calliptera analyses (all_CalMas) I have included these in the examples subdirectory.

If using a different dataset then first use PLINK to generate these files:

E.g. The example Lake Masoko files were generated as follows:
`plink --vcf ~/rds/rds-rd109-durbin-group/projects/cichlid/massoko/all_CalMas/masked_snps/sex/all_CalMas_variants_genome_clean.vcf --make-bed --set-missing-var-ids @:# --out all_CalMas`

Or generally (input vcf can be gzipped or uncompressed):
`plink --vcf dataset_name.vcf.gz --make-bed --set-missing-var-ids @:# --out dataset_name`

### Formatting Input File 

The input file containing test results to be plotted must be tab-delimited text file where each row provides information for a separate site in the genome.
It must contain min. 3 columns containing the following: 
- Site chromosome number
- Site position
- Test value <br>
Which columns contain each of these variables should be indicated using arguments `-x` or `--poscol`, `-y` or `--chrcol` and `-t` or `--testcol` respectively. <br>
This format and the default selections for the chromosome and position columns corresponds to the standard output from GWAS analyses with GEMMA.

## Usage
___

You must provide the full path to the plot_zoom executable script to run it (i.e. don't try to put it in your path, it won't work)

Executing plot_zoom without any options will print a list of arguments & descriptors:

### Essential Arguments
| Argument      |    Type    | Description |
| :-------------- |:----------:| :----:|
| `-i` <br>`--input`  | File | Path to input results file. |
| `-b` <br>`--plinkbed` | File | Path to PLINK bed file for LD calculations (prefix only).  |
| `-a` <br>`--annot`| File      |    Path to GTF annotation file. |
| `-o` <br>`--outpref` | String |  Path prefix for output files (Plot and PLINK LD calculation files). | 
| `-t` <br>`--testcol`  | Int | Index of column containing test statistic to be plotted (1-indexed). |

### Optional Arguments 
| Argument      |    Type    | Default | Description | 
| :-------------- |:----------:| :----: | ----- |
| `-h` <br>`--header`  | 0 or 1 | 1 | Indicates whether the line of the input file `-i` is a header (1) or not (0). | 
| `-x` <br>`--chrcol`  | Int | 1 | Index of column containing site chromosome number (default corresponds to GEMMA output).| 
| `-y` <br>`--poscol`  | Int | 3 | Index of column containing site position (default corresponds to GEMMA output). | 
|`-c` <br>`--fschr`  | Int | Chromosome number of first entry in input file | Focal SNP: **chromosome number** (default is top SNP only if input is sorted by statistic value). | 
|`-p` <br>`--fspos`  | Int | Position of first entry in input file | Focal SNP: **position**. Can also be used to specify other SNPs to be labelled in the plot by listing them, comma separated. (default is top SNP only if input is sorted by statistic value) | 
| `-w` <br> `--window` | Int,Int | 50000,50000 | Start and end positions defining the window to be plotted (default uses a window 50,000bp upstream and downstream of the focal SNP. | 
| `-s` <br> `--sigline` | Numeric or 'bfc' | NA |Value of significance threshold value to be plotted. Option bfc will plot a Bonferroni correction for the number of tests for p=0.05. If unused no sig line plotted. |
| `-l` <br> `--logtrans` | 0 or 1 | 1 |           Plot given statistic values (0) or their -log10 transform (1). <br>(default is to log transform)
| `-d` <br> `--diameter` | Float | 1 | Point diamater scaling factor. |
| `-f` <br> `--fontsize` | Float | 1 | Base font size scaling factor. | 
| `-g` <br> `--genelabs` | 0 or 1 | 1 | Include (1) or exclude (1) gene name annotations. |
| `-z` <br> `--highlight` | Int,Int | NA | Start and end genome positions of a region to highlight. |

### Ouput files

- Plot consisting of 3 panels (outpref_regional_plot.png) : SNP density plot; SNP test values coloured by LD with focal SNP; Gene Annoations
- Table of test values and R2 values (wrt focal SNP) for each plotted SNP in the window (outpref_plotted_SNP_values.txt)


### Example

`./plot_zoom.sh \` <br>
` -i masoko_sex_GWAS.assoc.txt \` <br>
` -t 14 \` <br>
` -b all_CalMas.bed` <br>
` -a Astatotilapia_calliptera.fAstCal1.2.99.chr.gtf \` <br>
` -c 19 \` <br>
` -p 21581905,21728957 \` <br>
` -w 21500000,21750000 \` <br>
` -s bfc \` <br>
` -o masoko_sex_GWAS_chr19`
#!/bin/bash 

# plot_zoom.sh

# Help Information

if [ "$#" -lt 1 ]
then
echo "--------------------------plot_zoom--------------------------"
echo "Description: Wrapper around PLINK and region_plotter.R to output local plots of the results of genomic analyses."
echo "Dependencies: PLINK, R: Bioconductor, rtracklayer"
echo "Usage: bash plot_zoom.sh [options]"
echo "ARGUMENTS"
echo "Required:"
echo "-i | --input      <string>        Path to input file containing results to be plotted."
echo "-t | --testcol    <int>           Number of column (1-indexed) containing test statistic to be plotted."
echo "-b | --plinkbed   <string>        Path to PLINK BED file (generated directly from VCF) to be used for LD calculations (prefix only)."
echo "-a | --annot      <string>        Path to GTF annotation file."
echo "-o | --outpref    <string>        Path/filename for output files (Plot and PLINK LD calculation files)."
echo "Optional:"
echo "-h | --header     <0|1>           First line of input file is not (0) or is (1) a header. (default assumes header present)"
echo "-x | --chrcol     <int>           Number of column (1-indexed) indicating site chromosome. (default is col 1)"
echo "-y | --poscol     <int>           Number of column (1-indexed) indicating site position. (default is col 3)"
echo "-c | --fschr      <int>           Focal chromosome. (default is firt entry in input file i.e. top SNP only if sorted by pval)"
echo "-p | --fspos      <int|int,-,->   Focal SNP position and pos of any other SNPs to be labelled, comma separated. (default focal SNP as above)"
echo "-w | --window     <int,int>       Plotting window start and end positions, comma separated. (default plots 50,000 bp either side of focal SNP)"
echo "-s | --sigline    <float|string>  Significance threshold value OR 'bfc' to plot Bonferroni corrected p=0.05. If unused no sig line plotted."
echo "-l | --logtrans   <0|1>           Plot absolute statistic values (0) or their -log10 transform (1). (default is to log transform)"
echo "-d | --diameter   <float>         Diameter of points (scaling factor) in plot. (default is 1)"
echo "-f | --fontsize   <float>         Base font size scaling factor. (default is 1)"
echo "-g | --genelabs   <0|1>           Setting to 0 removes gene annotation labels (default is to include them, 1)"
echo "-z | --highlight  <int,int>       Start and end positions of a region to highlight, command separated. (no highlight if option not supplied)"
exit
fi

### parse user arguments ###

while [[ $# -gt 1 ]]
	do
	key="$1"

case $key in
	-i|--input)
	IPUT="$2"
	shift
	;;
	-h|--header)
	HEAD="$2"
	shift
	;;
	-t|--testcol)
	TEST="$2"
	shift
	;;
	-x|--chrcol)
	CHRCOL="$2"
	shift
	;;
	-y|--poscol)
	PSCOL="$2"
	shift
	;;
	-a|--annot)
	ANNOT="$2"
	shift
	;;
	-b|--plinkbed)
	PBFILE="$2"
	shift
	;;
	-c|--fschr)
	CHROM="$2"
	shift
	;;
	-p|--fspos)
	POSLIST="$2"
	shift
	;;
	-w|--window)
	WIN="$2"
	shift
	;;
	-s|--sigline)
	SIG="$2"
	shift
	;;
	-l|--logtrans)
	LOG="$2"
	shift
	;;
	-d|--diameter)
	DIAM="$2"
	shift
	;;
	-f|--fontsize)
	FONT="$2"
	shift
	;;
	-g|--genelabs)
	GLAB="$2"
	shift
	;;
	-z|--highlight)
	HIL="$2"
	shift
	;;
	-o|--outpref)
	OUTPREF="$2"
	shift
	;;
	*)
	>&2 echo "Unknown command $1"
	exit 1
	;;
esac
shift
done

# set directory of plot_zoom script
my_dir="$(dirname "$0")"

## Print options being used

echo "--------------------------plot_zoom--------------------------"
echo "Description: Wrapper around PLINK and region_plotter.R to output regional plots of the results of genomic analyses."
echo "Dependencies: PLINK, R: Bioconductor, rtracklayer"
echo "Usage: bash plot_zoom.sh [options]"
echo "---------------------------OPTIONS----------------------------"

# CHECK INPUT FILE 
if [ -z "$IPUT" ]
then
     echo "** ERROR: No input results file given to be plotted. Option: -i | --input. **"
     exit
elif ! [ -f "$IPUT" ]
then
	 echo "** ERROR: Input results file can't be found. Option: -i | --input. **"
	 exit
else
     echo "INPUT: $IPUT"
fi

# CHECK PLINK BED FILE 
if [ -z "$PBFILE" ]
then
      PBFILE="BAD"
      echo "** ERROR: No PLINK bed file given to calculate LD. Option: -b | --plinkbed. **"
elif ! [ -f "${PBFILE}.bed" ]
then
	  PBFILE="BAD"
	  echo "** ERROR: PLINK bed file can't be found - please supply prefix ONLY (myfile not myfile.bed) . Option: -b | --plinkbed. **"
else
      echo "PLINKBED: $PBFILE"
fi

# CHECK ANNOTATION FILE
if [ -z "$ANNOT" ]
then
     ANNOT="BAD"
     echo "** ERROR: No GTF annotation file given. Option: -a | --annot. **"
elif ! [ -f "$ANNOT" ]
then
	 ANNOT="BAD"
	 echo "** ERROR: GTF annotation file can't be found. Option: -a | --annot. **"
else
     echo "ANNOTATION: $ANNOT"
fi

# CHECK HEADER ARGUMENT / SET DEFAULT 
if [ "$HEAD" == 0 ]
then
	 echo "HEADER: No."
elif [ "$HEAD" == 1 ]
then
	 echo "HEADER: Yes."
elif [ -z "$HEAD" ]
then
     HEAD=1
     echo "HEADER: Yes (default). Option -h | --header"
else
	HEAD="BAD"
	echo "** ERROR: Input given to option -h | --header is not 0 or 1. **"
fi

# CHECK TEST COLUMN IS VALID
# check number of columns in input file for comparison
ncol=$(awk -F'\t' '{print NF; exit}' "$IPUT")
if [ -z "$TEST" ]
then
     TEST="BAD"
     echo "** ERROR: No test statistic column specified for plotting. Option: -t | --testcol **" 
elif ! [[ "$TEST" =~ ^[1-9]+[0-9]*$ ]]
then
	 TEST="BAD"
	 echo "** ERROR: Input provided to test column is not a positive integer. Option -t | --testcol **"
elif [ "$TEST" -gt "$ncol" ]
then
     TEST="BAD"
     echo "** ERROR: Column $TEST selected to be plotted but only $ncol columns in input file. Option: -t | --testcol **"
else
	 echo "TEST: Plotting test statistic values from column $TEST"
fi

# CHECK CHROMOSOME COLUMN IS VALID
if [ -z "$CHRCOL" ]
then
	CHRCOL=1
     	echo "CHRCOL: 1 (default). Option: -x | --chrcol" 
elif ! [[ "$CHRCOL" =~ ^[1-9]+[0-9]*$ ]]
then
	  CHRCOL="BAD"
	 echo "** ERROR: Input provided to position column is not a positive integer. Option -x | --chrcol **"
elif [ "$CHRCOL" -gt "$ncol" ]
then
      CHRCOL="BAD"
      echo "** ERROR: Column $CHRCOL indicated as chromosome column but only $ncol columns in input file. Option: -x | --chrcol **"
else
	 echo "CHRCOL: $CHRCOL"
fi

# CHECK POSITION COLUMN IS VALID 
if [ -z "$PSCOL" ]
then
      PSCOL=3
      echo "POSCOL: 3 (default). Option: -y | --poscol" 
elif ! [[ "$PSCOL" =~ ^[1-9]+[0-9]*$ ]]
then
	 PSCOL="BAD"
	 echo "** ERROR: Input provided to position column is not a positive integer. Option -y | --poscol **"
elif [ "$TEST" -gt "$ncol" ]
then
	 PSCOL="BAD"
     echo "** ERROR: Column $PSCOL indicated as chromosome column but only $ncol columns in input file. Option: -y | --poscol **"
else
	 echo "POSCOL: $PSCOL"
fi

# CHECK FOCAL SNP 
# split 
IFS="," read -ra POSITION_LIST <<< "$POSLIST"
POS="${POSITION_LIST[0]}"

# default focal SNP

if [ -z "$CHROM" ] || [ -z "$POS" ]
then
	  CHROM="$(awk -vchr="$CHRCOL" 'NR==2 {print $chr}' $IPUT)"
	  POS="$(awk -vpos="$PSCOL" 'NR==2 {print $pos}' $IPUT)"
      echo "FOCALSNP: WARNING - Missing or incomplete focal SNP. Supply chromosome with option: -c | -fschr & position with option: -p | -fspos." 
      echo "          As a default the first SNP in the input file will be used - chr $CHROM pos $POS. (top SNP if input has been sorted by pval)"
else
      echo "FOCALSNP: chromosome $CHROM position $POS" 
fi 


# CHECK WINDOW 
# split input into separate variables 

IFS="," read -ra LIMS <<< "$WIN"
LEFT="${LIMS[0]}"
RIGHT="${LIMS[1]}"


if [ -z "$LEFT" ]
then
      LEFT="$(( ${POS} - 50000 ))"
elif ! [[ "$LEFT" =~ ^[0-9]+$ ]]
then 
	  LEFT="BAD"
fi 

if [ -z "$RIGHT" ]
then
      RIGHT="$(( ${POS} + 50000 ))"
elif ! [[ "$RIGHT" =~ ^[0-9]+$ ]]
then 
	  RIGHT="BAD"
fi 

if [ "$RIGHT" == BAD ] || [ "$LEFT" == BAD ]
then
      echo "** ERROR : Plotting window input ill-formatted. Option: -w | --window. **"
elif [ "$POS" -gt "$RIGHT" ] || [ "$POS" -lt "$LEFT" ] 
then
      echo "** ERROR : Focal SNP $POS is not within selected window of ${LEFT} - ${RIGHT}. **"
      exit
else
	  echo "PLOTTING WINDOW: chr${CHROM}:${LEFT}-${RIGHT}"
fi


# CHECK SIGNIFICANCE THRESHOLD
# perform test to see if significance threshold is interpretable as numeric by R, sig_isgood = 1 if good, 0 if bad. 
sig_isgood=$( echo $SIG | R --vanilla --slave -e 'rv=tryCatch(scan("stdin", quiet=TRUE), error=function(e){}); if (is.null(rv)) cat(0) else cat(1)')

if [ -z "$SIG" ]
then
      SIG="NA"
      echo "Sig: No significance threshold will be plotted. Option: -s | --sigline"
elif [ "$SIG" == "bfc" ]
then
	# Calculate bfc significance threshold in R 
	if [ "$SIG" == "bfc" ] && [ "$HEAD" -eq 1 ]
	then
		LINES=$(tail -n +2 $IPUT | wc -l)
		SIG=$( echo $LINES | R --vanilla --slave -e 'sites=as.numeric(scan("stdin", quiet=TRUE)); bfc=(0.05/sites); cat(bfc)')
	else
		LINES=$(wc -l $IPUT)
		SIG=$( echo $LINES | R --vanilla --slave -e 'sites=as.numeric(scan("stdin", quiet=TRUE)); bfc=(0.05/sites); cat(bfc)')
	fi 
	echo "SIG: A significance threshold of p=0.05 bonferroni corrected for total number of tests, calculated as ${SIG}, will be plotted."
elif [ "$sig_isgood" == 1 ]
then
	 echo "SIG: Significance threshold set at $SIG"  
else
	 SIG="BAD"
	 echo "** ERROR: Significance threshold supplied is not a valid numerical input. Option: -s | --sigline **" 
fi

# CHECK LOG TRANSFORMATION ARGUMENT / SET DEFAULT
if [ "$LOG" == 0 ]
then
	 echo "LOG: Plotting points WITHOUT log transformation."
elif [ "$LOG" == 1 ]
then
	 echo "LOG: Plotting -log10 transformation of points."
elif [ -z "$LOG" ]
then
     LOG=1
     echo "LOG: Plotting -log10 transformation of points by default. Option -l | --logtrans"
else
	 LOG="BAD"
	 echo "** ERROR: Input given to option -l | --logtrans is not 0 or 1. **"
fi


# CHECK DIAMETER ARGUMENT / SET DEFAULT 
if [ -z "$DIAM" ]
then
	DIAM=1
elif ! [[ "$DIAM" =~ ^([0-9]*[.])?[0-9]+$ ]]
then
	DIAM="BAD"
	echo "** ERROR: Input supplied to option -d | --diameter is not a float. **"
else
	echo "DIAMETER: Diamater of points in plot scaled by ${DIAM}"
fi

# CHECK FONTSIZE ARGUMENT / SET DEFAULT
if [ -z "$FONT" ]
then
	FONT=1
elif ! [[ "$FONT" =~ ^([0-9]*[.])?[0-9]+$ ]]
then
	FONT="BAD"
	echo "** ERROR: Input supplied to option -f | --fontsize is not a float. **"
else
	echo "FONTSIZE: Font sizes scaled by ${FONT}"
fi

# CHECK GENE LABELS ARGUMENT / SET DEFAULT 
if [ "$GLAB" == 0 ]
then
	 echo "LABEL GENES: No."
elif [ "$GLAB" == 1 ]
then
	 echo "LABEL GENES: Yes."
elif [ -z "$GLAB" ]
then
     GLAB=1
     echo "LABEL GENES: Yes (default). Option -g | --genelabs"
else
	GLAB="BAD"
	echo "** ERROR: Input given to option -g | --genelabs is not 0 or 1. **"
fi

# CHECK HIGHLIGHT / SET DEFAULT
if [ -z "$HIL" ]
then
	HIL="NA,NA"
else
	IFS="," read -ra HLIMS <<< "$HIL"
	LEFT_H="${HLIMS[0]}"
	RIGHT_H="${HLIMS[1]}"
	if ! [[ "$LEFT_H" =~ ^[0-9]+$ ]]
	then 
	  HIL="BAD"
	  echo "** ERROR : Highlight input ill-formatted. Option: -z | --highlight. **"
	elif ! [[ "$RIGHT_H" =~ ^[0-9]+$ ]]
	then 
	  HIL="BAD"
	  echo "** ERROR : Highlight input ill-formatted. Option: -z | --highlight. **"
	else
	  echo "HIGHLIGHT REGION: chr${CHROM}:${LEFT_H}-${RIGHT_H}"
	fi 
fi


# CHECK OUTPUT LOCATION PREFIX / SET DEFAULT
if [ -z "$OUTPREF" ]
then
	 OUTPREF="chr${CHROM}_${POS}_"
	 echo "OUTPUT: No output prefix given. Outputting files to working directory using default naming."
else
	 echo "OUTPUT: $OUTPREF"
fi

if [ "$PBFILE" == "BAD" ] || [ "$ANNOT" == "BAD" ] || [ "$TEST" == "BAD" ] || [ "$SIG" == "BAD" ] || [ "$LOG" == "BAD" ] || [ "$HEAD" == "BAD" ] || [ "$LEFT" == "BAD" ] || [ "$RIGHT" == "BAD" ] || [ "$DIAM" == "BAD" ] || [ "$FONT" == "BAD" ] || [ "$GLAB" == "BAD" ] || [ "$HIL" == "BAD" ] 
then
	 exit
else
	 echo "Ready to run."     
fi

# Convert window to spacing from focal SNP in kbfor PLINK
LEFT_SIZE=$(( $POS - $LEFT ))
RIGHT_SIZE=$(( $RIGHT - $POS ))
if [ "$LEFT_SIZE" -gt "$RIGHT_SIZE" ]
then
	WIDTH="$LEFT_SIZE"
else
	WIDTH="$RIGHT_SIZE"
fi

window_plink=$( echo $WIDTH | R --vanilla --slave -e 'width_bp=as.numeric(scan("stdin", quiet=TRUE)); width_kb=ceiling(width_bp/1000); cat(width_kb)')

# Generate temporary input file with only entries from chromosome being plotted 
awk -F, -v X="$CHROM" -v Y="$CHRCOL" '{if(int($Y) == X) {print}}' $IPUT > input_temp.txt
sed -i "1s/^/$(head -n1 $IPUT)\n/" input_temp.txt

# Generate temporary annotation file with only entries from chromosome being plotted
awk -F, -v X="$CHROM" '{if(int($1) == X) {print}}' $ANNOT > annot_temp.gtf

echo "---------------------PLINK_LD_CALCULATION---------------------"
# Calculate LD between focal SNP and surrounding SNPs
plink --bfile $PBFILE \
--r2 \
--ld-snp chr$CHROM:$POS \
--ld-window-kb $window_plink \
--ld-window 99999 \
--ld-window-r2 0 \
--out ${OUTPREF}

# Run plotting R script
Rscript $my_dir/region_plotter.R input_temp.txt "${OUTPREF}.ld" "$TEST" "$CHRCOL" "$PSCOL" "$CHROM" "$POSLIST" "$SIG" "$LOG" "$HEAD" "$OUTPREF" annot_temp.gtf "$LEFT" "$RIGHT" "$DIAM" "$FONT" "$GLAB" "$HIL"

# Remove temporary files
rm input_temp.txt
rm annot_temp.gtf
rm ${OUTPREF}.ld
rm ${OUTPREF}.log
rm ${OUTPREF}.nosex

exit
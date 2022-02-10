#@ Script for extracting lists of integration sites (unique read 1) and sisters (unique paires of read1-read2) from aligned data. 

# set up required paths
BAMDIR=<path/to/data>"/04_align"
OUTDIR=<path/to/data>"/05_extract_sites"
SCRIPTPATH=<path/to/scripts/dir>

# choose file to process from job variable
FILES=($(ls $BAMDIR/*/*/*.bam))
f=${FILES[0]}   # it is recommended to do in parallel (eg using batch job in hpc). 

echo $f

fbname=$(basename "$f" .bam)
dbname=$(dirname "$f")

# samtools view is used to filter the reads by mapq and exclude secondary mapping and upstream LTR reads. 
samtools view -q 10 -f 66 -F 256 $f | grep -v upstream > $dbname/$fbname"_R1.sam"
samtools view -q 10 -f 130 -F 256 $f | grep -v upstream > $dbname/$fbname"_R2.sam"

# if the file is empty, skip. Otherwise, process using R. 
lines1=`wc -l $dbname/$fbname"_R1.sam" | cut -d " " -f 1`
lines2=`wc -l $dbname/$fbname"_R2.sam" | cut -d " " -f 1`

if [ $lines1 -eq 0 ] || [ $lines2 -eq 0 ]
then
	echo -e "\n---\n"$fbname" - no reads remain\n----\n"
else
	Rscript --vanilla $SCRIPTPATH"/04_extractsites.R" $f
fi



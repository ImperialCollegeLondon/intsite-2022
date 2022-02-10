#@ This script will loop all files, and align them against a combined reference of hg38 and virus sequence.  

REFPATH=<path/to/reference/files>
FASTQPATH=<path/to/trimmed/files>
BWA=<path/to/bwa>
OUTPUTDIR=<path/to/aligned/output>


# choose file to process from job variable
FILES=($(ls $FASTQPATH/*val_1.fq.gz))
f=${FILES[0]}  # it is recommended to do in parallel (eg using batch job in hpc). 

fbname=$(basename "$f" _R1_val_1.fq.gz)
dbname=$(dirname "$f")

filename1=$dbname/$fbname"_R1_val_1.fq.gz"
filename2=$dbname/$fbname"_R2_val_2.fq.gz"

echo -e "$fbname\n----\n"


#### ALIGNMENT

refdir=<path/to/refernce/files/and/bwa/index>
outputfile=<path/to/output/bam>

# align using BWA
COMMAND="$BWA mem -M -t 8 \
    $REFPATH/$refdir \
    $filename1 $filename2 | samtools view -Shb - > $outputfile" 

echo -e "\n----\n$COMMAND\n----\n"

eval $COMMAND




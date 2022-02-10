#@ example script for preparing a fastq file. 
#@ This script may need to be prepared differently, depending on when a sequencing experiment was carried out, as different Illumina software version output the data differently (e.g. are the files pre demultiplexed)
#@ This script may need to be prepared differently, depending on the integration site sequencing approach (e.g. are the last 5 bases of the LTR reported)
#@ This script uses cutadapt to filter reads only to those which contain the last 5 bases of the HTLV/HIV LTR at the start of a separate read.  
#@ desired output - Fastq file, different file per barcode/lane combination, filtered for the last 5 bases of the LTR. 

RAWPATH=<path/to/raw/fastq/dir/for/single/flowcell>
OUTPUTPATH=<path/to/fastq/dir/for/single/flowcell>
CUTADAPT=<path/to/cutadapt/bin> # here used cutadapt version 1.18

LANE=1   # it is recommended to do in parallel (eg using batch job in hpc). 
FILES=$(echo -e $RAWPATH"/Unaligned_less5bases/*/*_L00"$LANE"_R1*fastq*")

for f in $FILES
do
    f2=`echo $f | sed 's#_R1_#_R2_#g' -`
    ffilter=`echo $f | sed 's#.fastq.gz#.fastq#g' - | sed 's#Unaligned_less5bases#Unaligned_5bases#g' -`

    fbname=$(basename "$f" _001.fastq.gz)

    vuindex=`echo $fbname |  awk '{n=split($0,a,"_");print a[n-2]}'`
    lane=`echo $fbname | awk '{n=split($0,a,"_");print a[n-1]}'`
        
    echo -e `date` $fbname

    output1=$OUTPUTPATH/{name}/$FLOWCELL"_"$vuindex"_"$lane"_R1.fq"
    output2=$OUTPUTPATH/{name}/$FLOWCELL"_"$vuindex"_"$lane"_R2.fq"

    COMMAND1="$CUTADAPT \
     -g HTLV1=^ACACA \
     -g HIV1=^TAGCA \
     -e 0 \
     --overlap 5 \
     --no-trim \
     --discard-untrimmed \
     --report=minimal \
     -o /homes/anatm/lmpcr/test/{name} \
     -p $output1 \
     $ffilter $f
    "

    COMMAND2="$CUTADAPT \
     -g HTLV1=^ACACA \
     -g HIV1=^TAGCA \
     -e 0 \
     --overlap 5 \
     --no-trim \
     --discard-untrimmed \
     --report=minimal \
     -o /homes/anatm/lmpcr/test/{name} \
     -p $output2 \
     $ffilter $f2
    "

    eval $COMMAND1
    eval $COMMAND2

done

echo -e `date` done!

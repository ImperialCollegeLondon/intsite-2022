#@ This script will loop across all files in a flowcell and use trimgalore to trim ilumina linkers and low quality reads.

FASTQPATH=<path/to/prepared/fastq> # see 01_prepare_fastq
OUTPUTPATH=<path/to/outputdir>

DIRS=($(ls -d $FASTQPATH/*/H*))
DIR=${DIRS[0]} # it is recommended to do in parallel (eg using batch job in hpc). 


# create output dir if does not yet exist. for example here based on flowcell code. 
fc=`echo $DIR | awk '{n=split($0,a,"/");print a[n-1]}'`
outputdir=`echo $OUTPUTPATH/$fc`

# report to user
echo 'Trim files with trim galore...'

# identify which sequence primer was used for trimming read2
# Note - trimming sequence is VIRUS SPECIFIC and will also differ if different library prep approach was used. 
read2trim=TGTGTACTAAGTT

echo read2 will be trimmed using: $read2trim



# trim files
# it is recommended to combine this step with FastQC or another application for QC before and after trim. 

for file in $DIR/*R1.fq
do
	echo $file

	if [ ! -s $file ]; then continue ; fi

	fbname=$(basename "$file" _R1.fq)
	read1=$fbname"_R1.fq"
	read2=$fbname"_R2.fq"

	# run trim galore on files
	trim_galore --paired \
	-q 20 \
	-o $outputdir \
	--stringency 3 \
	--gzip \
	-a AGATCGGAAGAGC -a2 $read2trim \
	--no_report_file $DIR/$read1 $DIR/$read2 &>> $outputdir/trimoutput.txt

done


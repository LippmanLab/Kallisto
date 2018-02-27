#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=4G
#$ -pe threads 8
#$ -l tmp_free=50G

#$ -o /sonas-hs/lippman/hpc/data/Pericarp_RNAseq/log
#$ -t 1-10

###
# Produced by Zachary H. Lemmon on February 7, 2018. This script is invoked in parallel to download from SRA to the TMPDIR (saves hdd space on HPC)
###

if [[ -z $SGE_TASK_ID ]] ; then SGE_TASK_ID=1 ; fi

# Set some default values
Parameters=$(sed -n -e "$SGE_TASK_ID p" ParameterFile)
SampleName=$( echo "$Parameters" | awk '{print $1}' )
SRR=$( echo "$Parameters" | awk '{print $2}' )
kallistoidx=$( echo "$Parameters" | awk '{print $3}' )
trimmomatic=$( echo "$Parameters" | awk '{print $4}' )
bootstraps=$( echo "$Parameters" | awk '{print $5}' )

if [[ ! -d $trimmomatic || -z $SRR || ! -f $kallistoidx || -z $SampleName ]] ; then
	echo -ne "\nParameterFile seems incorrectly formated. Please double check!\n\t*ParameterFile is tab delimited, field 1 = SampleName, field 2 = SRR######, field 3 = /path/to/kallistoidx, field 4 = /path/to/trimmomatic.jar, field 5 = bootstrapnumber|BLANK\n\n"
	exit
fi

if [[ -z $bootstraps ]] ; then bootstraps=0 ; fi

if [[ -z $TMPDIR ]] ; then TMPDIR=$PWD ; fi

echo -ne "\nStarting kallisto run. options set are:\n"
echo -ne "\tSampleName : $SampleName\n"
echo -ne "\tSRRaccession : $SRR\n"
echo -ne "\tkallisto index : $kallistoidx\n"
echo -ne "\ttrimmomatic path : $trimmomatic\n"
echo -ne "\tbootstraps : $bootstraps\n"
date
echo -ne "\n"

echo -ne "#############################\n#############################\n\n"
echo -ne "\nSRA accession provided ("$SRR"). \n\tChecking SRR syntax..." ;
if [[ "$SRR" =~ ^[ES]RR[0-9]{6,}$ ]] ; then
	echo -ne " correct!\n\n" ;
else
	echo -ne "incorrect SRA|ENA accession format, should be ^[ES]RR[0-9]{6,}$. Please check.\n\n"
	exit 1
fi
firstsix=${SRR:0:6}
echo -ne "Starting wget for $SRR\n"; date
wget -c -q -O "$TMPDIR"/"$SRR".sra ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"$firstsix"/"$SRR"/"$SRR".sra
echo -ne "wget complete\n" ; date
## Dump reads to split fastq.gz in the TMPDIR.
echo "running fastq-dump for $SRR" ; date
fastq-dump -O "$TMPDIR" --split-files --gzip "$TMPDIR"/"$SRR".sra
read1="$SRR"_1.fastq.gz
read2="$SRR"_2.fastq.gz
rm -f "$TMPDIR"/"$SRR".sra
echo -ne "fastq-dump completed for $SRR\n"; date

if [[ ! -f $TMPDIR/$read2 ]] ; then
	echo -ne "Read2 does not exist, indicating a single end run. Setting read2 to none\n"
	read2=none
fi

if [[ ! -f $TMPDIR/$read1 ]] ; then
	echo -ne "Expected read1 does not exist! Check fastq-dump\n"
	exit
fi

# Make some shorter names for the reads... cuts off everything before the last "/" of the path
read1short=${read1/*\/}
#echo $read1short
read2short=${read2/*\/}
#echo $read2short

# Quickly check read length... this is so we only keep reads that trim to at least half of the read length (awk '{print int($1/2)}')
if [[ $TMPDIR/$read1short =~ gz ]] ; then 
	echo "compressed file, running through gunzip -c to check length"
	minlen=$( gunzip -c $TMPDIR/$read1short | head -n2 | tail -n1 | wc -m | awk '{print int($1/2)}' )
else
	echo "uncompressed file, running through awk to check length..."
	minlen=$( head -n2 $TMPDIR/$read1short | tail -n1 | wc -m | awk '{print int($1/2)}' )
fi

# Run trimmomatic
if [ $read2 != "none" ] ; then
	echo "Trimming PE reads..." ; date
	java -jar "$trimmomatic/trimmomatic.jar" PE -threads 8 "$TMPDIR"/"$read1short" "$TMPDIR"/"$read2short" "$TMPDIR"/"$read1short"_P1.fastq "$TMPDIR"/"$read1short"_U1.fastq "$TMPDIR"/"$read2short"_P2.fastq "$TMPDIR"/"$read2short"_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15:1:FALSE LEADING:30 TRAILING:30 MINLEN:$minlen TOPHRED33
else
	echo "Trimming SE reads..." ; date
	java -jar "$trimmomatic/trimmomatic.jar" SE -threads 8 "$TMPDIR"/"$read1short" "$TMPDIR"/"$read1short"_P1.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:$minlen TOPHRED33
fi
echo "Done with trimmomatic..." ; date
	# Start kallisto run... pass through if else to run the correct PE or SE command.
if [ $read2 != "none" ] ; then
	echo "Starting PE kallisto run..." ; date
	kallisto quant --index=$kallistoidx --output-dir="$SampleName"_"$SRR" --threads=8 --bootstrap-samples=$bootstraps $TMPDIR/"$read1short"_P1.fastq $TMPDIR/"$read2short"_P2.fastq
else
	echo "Starting SE kallisto run..." ; date
	kallisto quant --index=$kallisto_idx --output-dir="$SampleName"_"$SRR" --single --fragment-length=200 --sd=40 --threads=8 --bootstrap-samples=$bootstraps $TMPDIR/"$read1short"_P1.fastq
fi
echo "Done with kallisto run..." ; date

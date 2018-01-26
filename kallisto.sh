#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=4G
#$ -pe threads 8
#$ -l tmp_free=30G

###
# Produced by Zachary H. Lemmon on November 3, 2017. This script is not done in parallel as it is very fast and thus workflow implemented copies samples from NLSAS one at a time so as not to overwhelm I/O capacity. Also doesn't require lots of free space on HPC for copying the raw data files over as they are directly copied to the $TMPDIR and processed there.
###

trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar # Path to trimmomatic jar
#kallisto_idx=$HOME/indexes/SL3.0/ITAG3.20/ITAG3.20_CDS_kallisto_index # Path to the kallisto index, make the index before run.
KALLISTOINDEX LINE
#samplelist=/sonas-hs/lippman/hpc/data/bop_RNAseq/Kallisto/Samples # pathtoread1	pathtoread2	samplename
SAMPLELIST LINE

while read read1 read2 sample rest; do
	echo -ne "\n#########################################\n#########################################\n\nStarting sample: $sample\n\n" ;

	# Make some shorter names for the reads... cuts off everything before the last "/" of the path
	read1short=${read1/*\/}
	#echo $read1short
	read2short=${read2/*\/}
	#echo $read2short

	# Copy the first read, and if read2 is not "none" then also copy that one.
	echo "Copying files over to temporary space." ; date
	cp -v $read1 $TMPDIR/$read1short
	#echo "willcp $read1 TMPDIR/$read1short"
	if [ $read2 != "none" ]	; then
		cp -v $read2 $TMPDIR/$read2short
		#echo "willcp $read2 TMPDIR/$read2short"
	else
		echo "This is a single end run. Not copying over second read."
	fi

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
		java -jar "$trimmomatic" PE -threads 8 "$TMPDIR"/"$read1short" "$TMPDIR"/"$read2short" "$TMPDIR"/"$read1short"_P1.fastq "$TMPDIR"/"$read1short"_U1.fastq "$TMPDIR"/"$read2short"_P2.fastq "$TMPDIR"/"$read2short"_U2.fastq ILLUMINACLIP:$HOME/bin/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:40:15:1:FALSE LEADING:30 TRAILING:30 MINLEN:$minlen TOPHRED33
	else
		echo "Trimming SE reads..." ; date
		java -jar "$trimmomatic" SE -threads 8 "$TMPDIR"/"$read1short" "$TMPDIR"/"$read1short"_P1.fastq ILLUMINACLIP:$HOME/bin/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:$minlen TOPHRED33
	fi
	echo "Done with trimmomatic..." ; date

	# Start kallisto run... pass through if else to run the correct PE or SE command.
	if [ $read2 != "none" ] ; then
		echo "Starting PE kallisto run..." ; date
		kallisto quant --index=$kallisto_idx --output-dir=$sample --threads=8 $TMPDIR/"$read1short"_P1.fastq $TMPDIR/"$read2short"_P2.fastq
	else
		echo "Starting SE kallisto run..." ; date
		kallisto quant --index=$kallisto_idx --output-dir=$sample --single --fragment-length=200 --sd=40 --threads=8 $TMPDIR/"$read1short"_P1.fastq
	fi
	echo "Done with kallisto run..." ; date

	#clear out temporary disk space for the next sample
	rm -fv $TMPDIR/$read1short $TMPDIR/$read2short $TMPDIR/"$read1short"_P1.fastq $TMPDIR/"$read2short"_P2.fastq $TMPDIR/"$read1short"_U1.fastq $TMPDIR/"$read2short"_U2.fastq

done < $samplelist

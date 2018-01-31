#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=4G
#$ -pe threads 8
#$ -l tmp_free=30G

###
# Produced by Zachary H. Lemmon on November 3, 2017. This script is not done in parallel as it is very fast and thus workflow implemented copies samples from NLSAS one at a time so as not to overwhelm I/O capacity. Also doesn't require lots of free space on HPC for copying the raw data files over as they are directly copied to the $TMPDIR and processed there.
###

# Set some default values
bootstraps=0

for i in "$@"; do
	case $i in
		--trimmomatic=*)
		trimmomatic="${i#*=}"
		shift
	;;
		--samplelist=*)
		samplelist="${i#*=}"
		shift
	;;
		--kallistoidx=*)
		kallistoidx="${i#*=}"
		shift
	;;
		--bootstraps=*)
		bootstraps="${i#*=}"
		shift
	;;
		--default)
		DEFAULT=YES
		shift
	;;
	esac
done

if [[ ! -f $trimmomatic || ! -f $samplelist || ! -f $kallistoidx ]] ; then
	echo -ne "\nOne or more of the options supplied don't lead to actual files! Check that all are correct\n\nusage: kallisto.sh --trimmomatic=/path/to/trimmomatic.jar --samplelist=/path/to/samplelist --kallistoidx=/path/to/kallistoidx/\n\t*samplelist file is tab delimited, field 1 = read1, field 2 = read2 (or NA), field 3 = SampleName\n\n"

echo "options set are:"
echo -ne "\ttrimmomatic path : $trimmomatic\n"
echo -ne "\tsamplelist : $samplelist\n"
echo -ne "\tkallisto index : $kallistoidx\n"
echo -ne "\tbootstraps : $bootstraps\n"

	exit
fi

if [[ -z $TMPDIR ]] ; then TMPDIR=$PWD ; fi

echo -ne "\nStarting kallisto run. options set are:\n"
echo -ne "\ttrimmomatic path : $trimmomatic\n"
echo -ne "\tsamplelist : $samplelist\n"
echo -ne "\tkallisto index : $kallistoidx\n"
echo -ne "\tbootstraps : $bootstraps\n"
date

while read read1 read2 sample rest; do
	echo -ne "\n#########################################\n#########################################\n\nStarting sample: $sample\n\n" ;

	if [[ ! -f $read1 || ! -f $read2 ]] ; then
		echo -ne "reads in samplelist do not exist! Check the path is correct\n"
		continue
	fi

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
		kallisto quant --index=$kallistoidx --output-dir=$sample --threads=8 $TMPDIR/"$read1short"_P1.fastq $TMPDIR/"$read2short"_P2.fastq -b=$bootstraps
	else
		echo "Starting SE kallisto run..." ; date
		kallisto quant --index=$kallisto_idx --output-dir=$sample --single --fragment-length=200 --sd=40 --threads=8 $TMPDIR/"$read1short"_P1.fastq -b=$bootstraps
	fi
	echo "Done with kallisto run..." ; date

	#clear out temporary disk space for the next sample
	rm -fv $TMPDIR/$read1short $TMPDIR/$read2short $TMPDIR/"$read1short"_P1.fastq $TMPDIR/"$read2short"_P2.fastq $TMPDIR/"$read1short"_U1.fastq $TMPDIR/"$read2short"_U2.fastq

done < $samplelist

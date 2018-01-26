#!/bin/bash

if [ $# != 3 ] ; then
	echo -e 'usage: ./prepKallisto.sh ProjectName /path/to/Samples /path/to/kallisto_index\n\tExample:\n\t\tProjectName = bopRNAseq\n\t\t/path/to/Samples = /sonas-hs/lippman/hpc/data/bop_RNAseq/Kallisto/Samples\n\t\t/path/to/kallisto_index = $HOME/indexes/SL3.0/ITAG3.20/ITAG3.20_CDS_kallisto_index'
	exit
fi

# Read in variables
ProjName=$1
path_to_samples=$2
path_to_kallistoidx=$3

# Modify kallisto.sh file
echo "Making kallisto_$ProjName.sh"
awk -v samples="$path_to_samples" -v kalidx="$path_to_kallistoidx" '$0 ~ /^SAMPLELIST LINE$/{print "samplelist=" samples; next}$0 ~ /^KALLISTOINDEX LINE$/{print "kallisto_idx=" kalidx; next}{print}' kallisto.sh > kallisto_$ProjName.sh
chmod a+x kallisto_$ProjName.sh

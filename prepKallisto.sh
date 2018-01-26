#!/bin/bash

if [ $# != 2 ] ; then
	echo -e 'usage: ./prepKallisto.sh /path/to/Samples /path/to/kallisto_index\n\tExample:\n\t\t/path/to/Samples = /sonas-hs/lippman/hpc/data/bop_RNAseq/Kallisto/Samples\n\t\t/path/to/kallisto_index = $HOME/indexes/SL3.0/ITAG3.20/ITAG3.20_CDS_kallisto_index'
	exit
fi

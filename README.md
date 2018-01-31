# Kallisto
## Nov 3, 2017
This is a pipeline implemented by Zachary H. Lemmon (<zlemmon@cshl.edu>). It uses the Kallisto pseudo-alignment algorithm to quantify expression of RNAseq samples. Runs samples in a series as Kallisto is a fast algorithm and you spend too much time copying files over to the HPC if you want to do in parallel. Doing samples in series also eliminates the requirement for large HPC free space to copy the raw fastq files over.

### Outline for running on a single node interactively:
1. Prepare 'Samples' file with three fields tab delimited in format: "FullpathtoRead1\tFullpathtoRead2\tSampleName\n"
2. kallisto.sh --trimmomatic=/path/to/trimmomatic.jar --samplelist=/path/to/samplelist --kallistoidx=/path/to/kallistoidx --bootstraps=NUMBERBOOT 

### Outline for running qsub:
1. Prepare 'Samples' file with three fields tab delimited in format: "FullpathtoRead1\tFullpathtoRead2\tSampleName\n"
2. qsub -pe threads 8 -l m_mem_free=4G -l tmp_free=40G -cwd -j y kallisto.sh

# Kallisto
## Nov 3, 2017
This is a pipeline implemented by Zachary H. Lemmon (<zlemmon@cshl.edu>). It uses the Kallisto pseudo-alignment algorithm to quantify expression of RNAseq samples. Runs samples in a series as Kallisto is a fast algorithm and you spend too much time copying files over to the HPC if you want to do in parallel. Doing samples in series also eliminates the requirement for large HPC free space to copy the raw fastq files over.

### Outline:
1. Prepare 'Samples' file with three fields tab delimited in format: "FullpathtoRead1\tFullpathtoRead2\tSampleName\n"
2. run './prepKallisto.sh /path/to/Samples /path/to/kallisto_index' 
3. 'qsub kallisto.sh' 

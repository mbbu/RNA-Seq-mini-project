#! usr/bin/bash
set -e

# Download and unzip gtf file
wget ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
gunzip Homo_sapiens.GRCh38.102.gtf.gz

# create a directory to save the counts
mkdir feature_Counts

# Generate counts	
featureCounts -a Homo_sapiens.GRCh38.102.gtf  -o feature_Counts/counts.txt ./hisat/*.bam


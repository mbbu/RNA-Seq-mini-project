#! /usr/bin/bash
set -e
mkdir RNA_Seq_Analysis
cd RNA_Seq_Analysis
#download the reads and metadata
sample_id=$'sample37\nsample38\nsample39\nsample40\nsample41\nsample42'
for sample in `echo "$sample_id"`
do
	wget -c http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${sample}_R1.fastq.gz
	wget -c http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/${sample}_R2.fastq.gz

done

wget http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/practice.dataset.metadata.tsv

#Quality check the raw reads using fastqc

# First create the directory where quality reports will be directed
mkdir quality_reports
mkdir quality_recheck_reports  #reports generated after sequence trimming will be saved in this directory

fastqc -o quality_reports *.fastq.gz #fastqc command

mkdir trimmomatic_result

# Trim reads with low quality and remove adapters if present using trimmomatic
 
for R1 in *R1* 
do
	R2=${R1//R1.fastq.gz/R2.fastq.gz}
	R1_paired=${R1//.fastq.gz/_paired.fastq.gz}
	R1_unpaired=${R1//.fastq.gz/_unpaired.fastq.gz}
	R2_paired=${R2//.fastq.gz/_paired.fastq.gz}
	R2_unpaired=${R2//.fastq.gz/_unpaired.fastq.gz}
	trimmomatic PE $R1 $R2 ./trimmomatic_result/$R1_paired ./trimmomatic_result/$R1_unpaired ./trimmomatic_result/$R2_paired ./trimmomatic_result/$R2_unpaired SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 HEADCROP:10 
done

# Perform quality recheck after trimming using fastqc

fastqc -o  ./quality_recheck_reports trimmomatic_result/*_paired.fastq.gz

# pseudo-alignment and generate gene counts using kallisto


time bash kallisto.sh

# Perform alignment using hisat2 and gerenate counts using feature counts

#time bash hisat2.sh

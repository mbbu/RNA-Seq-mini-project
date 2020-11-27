#! /usr/bin/bash

#Quality check the raw reads using fastqc

# First create the directory where quality reports will be directed
mkdir ../quality_reports
mkdir ../quality_reports/quality_recheck_reports/  #reports generated after sequence trimming will be saved in this directory

fastqc -o ../quality_reports *.fastq #fastqc command


# Trim reads with low quality and remove adapters if present using trimmomatic
 
for R1 in *R1* 
do
	R2=${R1//R1.fastq/R2.fastq}
	R1_paired=${R1//.fastq/_paired.fastq}
	R1_unpaired=${R1//.fastq/_unpaired.fastq}
	R2_paired=${R2//.fastq/_paired.fastq}
	R2_unpaired=${R2//.fastq/_unpaired.fastq}
	trimmomatic PE $R1 $R2 $R1_paired $R1_unpaired $R2_paired $R2_unpaired -summary summary.txt -trimlog logfile SLIDINGWINDOW:4:20 MINLEN:20 ILLUMINCACLIP:NexteraPE-PE.fa:2:40:15 
done

# Perform quality recheck after trimming using fastqc

fastqc -o  ../quality_reports/quality_recheck_reports *_paired.fastq

# perform pseudo-alignment and generate gene counts using kallisto

# create a kallisto index

kallisto index  -i kallisto_index Homo_sapiens.GRCh38.cdna.abinitio.fa.gz

#now run kallisto
for sample in *R1_paired.fastq*
do
	R2=${sample//_R1_paired.fastq/_R2_paired.fastq}
	kallisto quant -i ../Kallisto/kallisto_index -b 30 -o ./  ${sample} $R2
done



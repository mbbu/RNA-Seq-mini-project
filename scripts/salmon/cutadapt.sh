#! usr/bin/bash

set -e

#Trimming adapters from your sequences depending on the quality as assessed from fastqc.

mkdir cutadaptresults

for file in ./*_R2.fastq.gz
do
    
    file2=${file/_R2.fastq.gz/_R1.fastq.gz} # replace _R1.fastq.gz in the file name with _R2.fastq.gz
    out1=./cutadaptresults/${file}
    out2=${out1/_R2.fastq.gz/_R1.fastq.gz}
    cutadapt -m 20 -u 15 -q 30,30 -o ${out1} -p ${out2}  ${file} ${file2}
done

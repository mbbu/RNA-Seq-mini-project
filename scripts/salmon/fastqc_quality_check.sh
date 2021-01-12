#! usr/bin/bash

set e-

# create  a directory where your reports will be stored

mkdir qualityreports

# check the quality of the raw reads using fastqc

fastqc -o ./qualityreports *.fastq.gz




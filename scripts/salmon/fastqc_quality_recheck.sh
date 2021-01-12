#! usr/bin/bash

set -e

#After trimming do a quality re-check before proceeding to the next steps.

mkdir quality-re-check

#quality-re-check command

fastqc -o ./quality-re-check ./cutadaptresults/*.fastq.gz


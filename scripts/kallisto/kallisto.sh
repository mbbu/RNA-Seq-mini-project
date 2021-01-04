# create a kallisto index
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz


kallisto index gencode.v36.transcripts.fa.gz -i ./kallisto_index


#mkdir kallisto_results
#now run kallisto
for sample in ./trimmomatic_result/*R1_paired.fastq.gz*
do
	output=${sample//_R1_paired.fastq.gz/ }
	R2=${sample//_R1_paired.fastq.gz/_R2_paired.fastq.gz}
	kallisto quant -i kallisto_index -b 100 -o ./kallisto_results/$output ./trimmomatic_result/${sample} ./trimmomatic_result/$R2
done


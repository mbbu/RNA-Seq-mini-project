#! usr/bin/bash
set -e

#unzip the file to be indexed

wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


#Building index file

hisat2-build -f Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa.idx

#make a directory to save the index
mkdir -p hisat_index 
mv *.ht2 ./hisat_index
mkdir -p hisat
# Run hisat2 allignment
cd trimmomatic_result
for R1 in *R1_paired.fastq.gz*
do

	m2=${R1//R1_paired.fastq.gz/R2_paired.fastq.gz}
	hisat_output=${R1//R1_paired.fastq.gz/hisat.sam}

	hisat2 -x ../hisat_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa.idx -q -1 $R1 -2 $m2 -S ../hisat/$hisat_output


# convert to binary, sort and index

	samtools view -b ../hisat/$hisat_output | samtools sort  > ../hisat/${hisat_output//_hisat.sam/_hisat_sorted.bam} 
	samtools index ../hisat/${hisat_output//_hisat.sam/_hisat_sorted.bam} #./hisat/${hisat_output//_hisat.sam/_index}

#delete .sam
	rm $hisat_output

done
cd ../
bash featureCounts.sh

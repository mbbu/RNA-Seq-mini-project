#! usr/bin/bash

set -e

#Salmon quantification

#make directory where the quantifications will be stored

mkdir quants

#salmon  requires a decoy aware transcriptome index when performing quantification in mapping based mode . 
#Salmon indexing requires the names of the genome targets, which is extractable by using the grep command.

grep "^>" <(gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

#Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. 
#The genome targets (decoys) should come after the transcriptome targets in the reference

cat gencode.v36.transcripts.fa.gz Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > congentra.fa.gz

#Building the index

salmon index -t congentra.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode #--gencode flag is for removing extra metadata in the target header separated by | from the gencode reference. 

# salmon loop

for file in *_R2.fastq.gz
do  
    samp=${file//_R2.fastq.gz/_R1.fastq.gz}
    echo "Processing ${file}"
    salmon quant -i salmon_index -l A \
         -1 ${file}\
         -2 ${samp}\
         -p 2 --validateMappings -o ./quants/${file//_R2.fastq.gz/} 
done





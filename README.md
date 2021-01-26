# RNA-Seq-mini-project

RNA-seq (RNA-sequencing) is a technique that can examine the quantity and sequences of RNA in a sample using next-generation sequencing (NGS). Over the past few years, RNA sequencing (RNA-seq) has become an indispensable tool for transcriptome-wide analysis of differential gene expression and differential splicing of mRNAs. It is rapidly replacing gene expression microarrays in many labs as it lets you quantify, discover, and profile RNAs.

Several tools and pipelines exist for RNA-Seq data analysis. Different consortiums and institutions use different sets of guidelines and standards for their data analysis. The [H3ABioNet](https://www.h3abionet.org) has developed a standard SOP and guidelines for RNA-Seq data analysis with some recommendations for gene expression analysis in human.

In this repo, we document RNA-seq data analysis following this [guidelines](https://h3abionet.github.io/H3ABionet-SOPs/RNA-Seq) developed by H3ABioNet. Data used in this project is available [here](http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/)

## Roadmap

**Phase I (Pre-processing analysis),  Time: 1 week.**

 *Tools : fastqc - 01109 , Trimmomatic - 0.39v*

  - Download raw reads
  - Check quality of the raw reads
  - Adapter removal and quality trimming
  - Quality recheck

**Phase II (Gene Expression Analysis ),  Time: 2 Days.**

 Generate gene/transcript level counts 

 *Tool - kallisto*
  
 - Align reads to reference genome 
 - Generate estimated counts using pseudo-alignment approach
 - Collecting and tabulating alignment stats

**Phase III (R - Analysis ),  Time: 4 Days**

 - QC and outlier removal / Batch detection.
 - Answer general questions of the project
 - wrap-up


## Setup

Create conda environment

    $ conda create --name [environ_name]
    
Activate conda environment

    $ conda activate [environ_name]
    
Install tools

    $ conda install [toolname] -c bioconda

## Tools:

| Tool name    |  Version        |   Use       |
|------|-----------|------------------------|
|  [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)            |      0.11.9           |    Check the quality of the reads          |
|  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)            |    0.39             |    Trim adapter remnants and low quality reads          |
| [Kallisto]( https://pachterlab.github.io/kallisto/)                   |  0.46.2                     |      pseudo-alignment and gene counts                      

**R-Analysis**                            
 |  Package           | Use|
  |--------------------|----|
  |DESeq2| To analyse count data and test for differential expression.|
|rhdf5 | To read abundance.h5 file |
|tximport| To import abundance.h5 file|
|pheatmap|  To draw clustered heatmaps |
|RcolorBrewer| Contains a ready-to-use color palettes for creating heatmaps|
|tximportData| Provides output of running Kallisto|





## Workflow:
**Phase 1**
- [x] Download raw reads
- [x] Quality check of the raw reads
- [x] Adapter removal and quality trimming
- [x] Quality recheck



**Phase 2**
- [x] Alignment
- [x] Trascripts/gene counts
- [x] Collect and tabulate statistics



**Phase 3**
- [x] Statistical analysis:
    - [x]  QC check
    - [x]  Outlier removal and normalization
    - [x]  Differential expression

**How to use the provided scripts for analysis**

**Hisat pipeline**

-The documents are found [here](https://github.com/mbbu/RNA-Seq-mini-project/tree/main/scripts/hisat2)

-First, put your raw reads and metadata in one file.In case of HPC make sure you ```module load``` all the tools required for this pipeline.
You will begin with checking the quality of your reads using [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).Here you will get the information on which reads to trim or not. Those that require trimming to remove low quality reads and reads that have a shorter length than your preffered length 
will proceed for trimming using  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). This was done using this [script](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/hisat2/fastqc-trimmomatic.sh)

-Allignment of the reads requires a reference genome that will used to create an index to be used for the allignment. Using the command ```wget``` you can obtain the fasta file in relation to your reads and use hisat2 in creation of indeces and proceed for allignment of the reads.This was done using this [script](https://github.com/mbbu/RNA-Seq-mini-project/tree/main/scripts/hisat2)

-When using HISAT2 the counts are obtained using a different tool.  In this pipeline, we used *features count* to count the reads that alligned to the indexes created from the reference genome.The counts were done using this [script](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/hisat2/featureCounts.sh)

-The counts obtained from featuresCount were used for statistical analysis in R using DESeq2. The statistical analysis done are contained in the [DESeq2 Rmd](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/hisat2/Features_R-analysis.Rmd).


**Salmon Pipeline Scripts**

The scripts are found [here](https://github.com/mbbu/RNA-Seq-mini-project/tree/main/scripts/salmon)

**Data**

Raw data/reads from the sequencer, Metadata(Downloaded from [here](http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/)), Reference genome,downloaded  
[here](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz) and the [Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](http://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) 
 in one directory as the scripts.


**Phase I (Pre-processing)**

Quality control check was done using this, [fastqc_quality_check.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/fastqc_quality_check.sh) script.

Data cleaning involves , removal of adapter remnants, short reads and low quality bases. Cutadapt trimming tool was preffered and the script [cutadapt.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/cutadapt.sh) was used.

Quality-recheck after trimming is necessary to examine the extent to which your data was cleaned and this was achieved using [fastqc_quality_recheck.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/fastqc_quality_recheck.sh) script.


**Phase II (Gene Expression Analysis)**

Involves Alignment of reads, Gene counts and Tabulating of the statistics, the script [salmon.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/salmon.sh) was used.

**Phase III (Statistical analysis/Differential Expression)**

EdgeR package was used for normalization, statistical analysis and visualization of the gene counts using this [EdgeR_Analysis_script.Rmd](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/EdgeR_Analysis_script.Rmd) script.
The generated html document for the EdgeR can be assesed [here](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/EdgeR%20Analysis_script.html).

**Kallisto Pipeline**

The scripts can be accessed from this repo or [here](https://github.com/mbbu/RNA-Seq-mini-project/tree/main/scripts/kallisto).

**Data**

The code for downloading the required data is included in the scripts, since the data is huge it takes time to download the data. Incase you've already downloaded your data as from above then you can hash out the download codes from the scripts or if you wish to obtain raw data separately then use the below links.

Raw data reads and metadata can be downloaded from [here](http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/) incase you didn't download them from above pipelines. Reference genome can also be downloaded from [here](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz).

**Phase I (Pre-processing)**

Quality control check was done using fastqc tool which informed the data cleaning parameters in the next step.
Data cleaning involves , removal of adapter remnants, short reads and low quality bases. Trimmomatic  trimming tool was preffered in this pipeline. 

Quality-recheck after trimming is necessary to examine the extent to which your data was cleaned.This was also done using Fastqc.

The combined script for Fastqc and Trimmomatic with details are in [fastqc-trimmomatic.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/kallisto/fastqc-trimmomatic.sh) script.


**Phase II (Gene Expression Analysis)**

Involves Alignment of reads, Gene counts and Tabulating of the statistics, the script [kallisto.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/kallisto/kallisto.sh) was used to achieve this.

**Phase III (Statistical analysis/Differential Expression)**

DESeq2 package was used for normalization, statistical analysis and visualization of the gene counts using this [kallisto_Deseq_analysis.Rmd](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/kallisto/kallisto_Deseq_analysis.Rmd) script.


The generated html document for the DESeq2 can be assesed [here](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/kallisto/kallisto_Deseq_analysis.html)




The final report analysis report is available [here](https://mbbu.github.io/RNA-Seq-mini-project/reports/Rnaseq--mini-project-report--1-.html)





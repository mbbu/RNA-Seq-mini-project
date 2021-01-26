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

*Salmon related scripts*

**How to use the provided scripts for analysis**

**Salmon pipeline scripts**

**Data**

Raw data/reads from the sequencer, Metadata(Downloaded from above), Reference genome,downloaded  
[here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz) and the 
[Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz](ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
 in one directory as the scripts.


**Phase I (Pre-processing)**

Quality control check was done using this, [fastqc_quality_check.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/fastqc_quality_check.sh) script.

Data cleaning involves , removal of adapter remnants, short reads and low quality bases. Cutadapt trimming tool was preffered and the script [cutadapt.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/cutadapt.sh) was used.

Quality-recheck after trimming is necessary to examine the extent to which your data was cleaned and this was achieved using [fastqc_quality_recheck.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/fastqc_quality_recheck.sh) script.


**Phase II (Gene Expression Analysis)**

Involves Alignment of reads, Gene counts and Tabulating of the statistics, the script [salmon.sh](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/salmon.sh) was used.

**Phase III (Statistical analysis/Differential Expression)**

EdgeR package was used for normalization, statistical analysis and visualization of the gene counts using this [Salmon-EdgeR_script.Rmd](https://github.com/mbbu/RNA-Seq-mini-project/blob/main/scripts/salmon/Salmon-EdgeR_script.Rmd) script.


The final report analysis report is available [here](https://mbbu.github.io/RNA-Seq-mini-project/reports/Rnaseq--mini-project-report--1-.html)

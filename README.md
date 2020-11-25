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
| [Kallisto]( https://pachterlab.github.io/kallisto/)                   |  0.46.2                     |      pseudo-alignment and gene counts                                                  |
 






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
    - [ ] QC check
    - [ ] Outlier removal and normalization
    - [ ] Differential expression


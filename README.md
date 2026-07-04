# Welcome to my Bulk RNA-seq Tutorial

*Work in progress*

In this tutorial, you will find a step-by-step protocol for achieving a reliable, simple, and reproducible bulk RNA-seq analysis. The tutorial covers the following workflow steps:

- **Data acquisition** from NCBI SRA using the SRA Toolkit
- **Folder structure** setup
- **Secondary analysis** using Bash and Nextflow scripting, including:
    - Preprocessing: QC, trimming, post-trimming QC, and duplicate marking
    - Alignment
    - Raw count table generation per sample
    - Visualisation (IGV)
- **Tertiary analysis** (RStudio):
    - Count matrix generation
    - DESeq2 transformation
    - PCA
    - Volcano and heatmap plots

The tutorial starts with three dataset samples, each representing one specific experimental group, from the study [DOI: 10.1038/s41467-018-07329-0](https://www.nature.com/articles/s41467-018-07329-0). Each dataset consists of paired-end FASTQ files (R1 and R2), which are used for the secondary analysis. For the tertiary analysis, a larger set of raw counts from more samples per experimental group from the same study is used to enable DESeq2 to normalise and calculate differential expression.

The purpose of this tutorial is educational and to promote reproducibility.

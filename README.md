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

## 🔬 Workflow Overview
Overview of the bulk RNAseq pipeline in this repository.

<p align="center">
  <img src="images/Gemini_Generated_Image_l5yy5vl5yy5vl5yy.png" 
       width="75%">
</p>

- *Image generated in collaboration with Gemini (Google AI) via iterative prompting.*

---

## Tutorial structure

### 1️⃣ Part I – Preparation & setup  

1. Folder structure  
2. Finding an open dataset project   
3. Downloading paired-end reads datasets by **SRA Toolkit**   
4. Reference genome indexes from HISAT2
5. Conda environments

➡️ **Start here:**  
👉 [Part I – Preparation & setup](README_Part1_setup_bulkrnaseq.md)

---

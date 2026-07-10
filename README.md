# Welcome to my Bulk RNA-seq Tutorial

![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-blue)
![Status](https://img.shields.io/badge/status-work_in_progress-orange)

*Work in progress*

In this tutorial, you will find a step-by-step protocol for achieving a reliable, simple, and reproducible bulk RNA-seq analysis, using **Bash**, **Nextflow**, and **R**. Designed for teaching purposes with real public data.  
  
The tutorial covers the following workflow steps:

- **Part I**: **Setup and data preparation** from NCBI SRA using the SRA Toolkit
- **Part II**: **Preprocessing & Secondary analysis** using Bash scripting and Nextflow workflows, including:
    - Preprocessing: QC, trimming, post-trimming QC
    - Alignment (`HISAT2`)
    - Flagging duplicates (`Picard`)
    - Determination of strandedness
    - Gene-level paired-end read quantification (raw counts) using (`featureCounts`)
    - Visualisation (IGV)
- **Part III**: **Tertiary analysis** using **R**:
    - Count matrix generation
    - DESeq2 transformation
    - PCA
    - Volcano and heatmap plots
- **Part IV**: Pseudoalignment (`salmon`)
    - Alignment vs pseudoalignment

The tutorial starts with two representative RNA-seq samples from the study [DOI: 10.1038/s41467-018-07329-0](https://www.nature.com/articles/s41467-018-07329-0). Each dataset consists of paired-end (R1 and R2) FASTQ files, which are used for the preprocessing and secondary analysis.  
The same RNA-seq workflow is implemented twice: first using Bash script to illustrate each analysis step, and later using Nextflow to demonstrate workflow automation and reproducibility.  
For the tertiary analysis, a larger set of raw counts from more samples per experimental group from the same study is used to enable DESeq2 to normalise and calculate differential expression (DE).

> [!IMPORTANT]  
> **DESeq2** requires biological replicates to estimate dispersion. Analyzing 1 sample per condition (n=1) is statistically invalid for differential expression. The two samples used in **Part I** and **Part II** are strictly for demonstrating the QC, trimming, alignment, and counting pipeline (due to file size limits). The actual DE analysis uses a separate, larger count matrix (included in the repo) with at least 3 biological replicates per group.

The purpose of this tutorial is educational and to promote reproducibility.

## 🔬 Workflow Overview
Overview of the bulk RNA-seq pipeline in this repository.

<p align="center">
  <img src="images/Gemini_Generated_Image_l5yy5vl5yy5vl5yy.png" 
       alt="Bulk RNA-seq Workflow" 
       style="max-width: 800px; width: 100%;">
</p>

- *Image generated in collaboration with Gemini (Google AI) via iterative prompting.*

---

## Tutorial structure

### 1️⃣ Part I – Setup and data preparation  

1. Folder structure  
2. Find & download paired-end RNA-seq datasets   
3. Downloading pre-built HISAT2 genome indexes (e.g., *Homo sapiens* GRCh38/hg38)   
5. Conda environments  

➡️ **Start here:**  
👉 [Part I – Setup & data preparation](README_Part1-3_setup_bulkrnaseq.md)

---

## 🧪 Tested Environment

- macOS (Intel)
- 8 GB RAM
- Free space: < 40GB
- macOS Big Sur 11.7.11 (Intel)
- Conda-based installation

> [!NOTE]  
> ≥ 16 GB RAM (or even 32 GB) would be recommended for a full-scale dataset analysis. An 8 GB RAM laptop may suffice only for small test datasets such as those used in this tutorial.


This repository is intended for **educational and research purposes only**.  
It is **not validated for clinical diagnostic use** and should not be used for medical decision-making.

All analyses are performed on publicly available research datasets.

---

## 🧬 Target audience

This tutorial is intended for:
- Bioinformatics students with or without biology or engineer backgrounds
- Researchers transitioning from wet lab to computational analysis
- Researchers who want a **transparent, Bash/Nextflow workflow** that can be run locally

Basic familiarity with Linux command line is recommended.

---

## 🔎 Third-Party Tools & Resources

This tutorial uses and displays output or screenshots generated from the following tools and databases:

- Integrative Genomics Viewer (IGV, Broad Institute)
- FastQC
- MultiQC
- NCBI Sequence Read Archive (SRA)

All trademarks, software, and database contents belong to their respective owners.  
Screenshots and outputs are shown for educational and demonstration purposes only.

---

## 📜 License

© 2026 **bioinfo-frano**

This project is licensed under the **MIT License**. See the full license [here](LICENSE).


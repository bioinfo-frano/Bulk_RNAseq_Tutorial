# Part II – Preprocessing & Secondary analysis

## Table of Contents

- [Introduction](#introduction)
- [Creating the computing environment for bulk RNA-seq analysis](#creating-the-computing-environment-for-bulk-rna-seq-analysis)  
    - [I. Create a folder structure](#i-create-a-folder-structure)  
    - [II. Find & download paired-end RNA-seq datasets](#ii-find--download-paired-end-rna-seq-datasets)  
    - [III. Download a pre-built HISAT2 genome index](#iii-download-a-pre-built-hisat2-genome-index)  
    - [IV. Create a Conda environment](#iv-create-a-conda-environment)  
    - [V. Create a BED12 file](#v-create-a-bed12-file)  
    - [VI. Final folder structure: before starting bulk RNA-seq analysis](#vi-final-folder-structure-before-starting-bulk-rna-seq-analysis)

## Introduction

Analysing transcriptomic datasets independently can be quite difficult if you don't have the proper guidance and, importantly, enough patience, time, and computational resources. At the end, you would have to send your datasets to an external bioinformatician or try other options all of which may involve additional costs. This is the logical solution when the statistics, tables and plots are urgently needed for the submission of scientific manuscripts or when preparing seminars.  
The purpose of this tutorial is to show that you can independently analyse your data using a personal computer (e.g., laptop) or workstation, which has limited computational resources.  
For the sake of learning, how to analyse your RNA-seq datasets, ideally, you should have some basic knowledge on command line, bash scripting, R programming, and python. However, if you don't have it, don't worry, **learn by doing!**

I would strongly suggest the following tutorials, so that you train yourself in these topics.
  
- [Learn Mac Terminal Basics - macmostvideo - YouTube](https://www.youtube.com/watch?v=ZkoEHvG3GI8)   
- [Command Line Basics for Beginners - Full Course - freeCodeCamp.org - YouTube](https://www.youtube.com/watch?v=mABpAI-pCw0)    
- [MASTERING Command Prompt Basics! | Tutorial (for Windows) - Skill Foundry - YouTube](https://www.youtube.com/watch?v=QBWX_4ho8D4)   
- [Bash Scripting Tutorial for Beginners  - TechWorld with Nana - YouTube](https://www.youtube.com/watch?v=PNhq_4d-5ek)   

Optional, but highly recommended too:   
  
- [R programming in one hour - a crash course for beginners - R Programming 101 - YouTube](https://www.youtube.com/watch?v=eR-XRSKsuR4&t=176s)   
- [Python Full Course for Beginners - Programming with Mosh - YouTube](https://www.youtube.com/watch?v=K5KVEU3aaeQ)   
  
> [!NOTE]  
> This guide was developed and tested on macOS running on Intel processors. Users on Apple Silicon (M1/M2/…/M5) or Linux systems may need to adapt certain steps.
  
<br> 

> [!IMPORTANT]
> **Conda prerequisites:**   
> This guide assumes that Miniconda3 is already installed on your computer. If not, please consult the official [documentation](https://docs.conda.io/projects/conda/en/stable/index.html) or watch this [YouTube](https://www.youtube.com/watch?v=hDGSZMLS5F4&t=67s) video.  
> When Miniconda is already installed, you should see the `(base)` environment activated in your Terminal.
  
<br>    



## VI. Final folder structure: before starting bulk RNA-seq analysis

```bash
Bulk_rnaseq/
├── data
│   └── PRJNA437330
│       ├── SRR6815993
│       │   └── raw_fastq
│       │       └── SRR6815993_{1,2}.fastq.gz
│       └── SRR6816017
│           └── raw_fastq
│               └── SRR6816017_{1,2}.fastq.gz
├── reference
│   ├── hisat2_index
│   │   ├── grch38_tran
│   │   │   └── genome_tran.{1,2,3,4,5,6,7,8}.ht2
│   │   └── grch38_tran.tar.gz
│   └── intervals
│       ├── gencode.v38.annotation.gtf.gz
│       ├── gencode.v38.annotation.bed
│       ├── gencode.v38.annotation.nochr.bed
│       └── gencode.v38.annotation.nochr.clean.bed
└── scripts
    └── RNA1_environment.yml
```


---

If you have reached the end of **PART I**, I congratulate you!!  
Continue to the 👉 [Part II – Secondary analysis](README_Part2-3_secondary_bulkrnaseq.md), where you'll start with the preprocessing analysis to alignment till the generation of raw counts tables, using bash and nextflow scripting explained step-by-step.

Go back to the top of 👉 [Part I – Setup & data preparation](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation)

Go to the main page 👉 [Bulk RNA-seq Tutorial](README.md)
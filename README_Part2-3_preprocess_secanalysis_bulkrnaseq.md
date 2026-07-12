# Part II – Preprocessing & Secondary analysis

## Table of Contents

- [Introduction](#introduction)  
- [Pipeline overview](#pipeline-overview)  
- [Bash pipeline](#bash-pipeline)
    - [I. Bash: Preprocessing](#i-bash-preprocessing)  
    - [II. Bash: Alignment and mark duplicates](#ii-bash-alignment-and-mark-duplicates)  
    - [III. Bash: Gene-level paired-end read quantification](#iii-bash-gene-level-paired-end-read-quantification)  

    - [I. Nextflow: Preprocessing](#i-nextflow-preprocessing)  
    - [II. Nextflow: Alignment and mark duplicates](#ii-nextflow-alignment-and-mark-duplicates)  
    - [III. Nextflow: Gene-level paired-end read quantification](#iii-nextflow-gene-level-paired-end-read-quantification)  
- [Visualization](#visualization)


## Introduction

In this second part of bulk RNA-seq analysis, the two datasets downloaded in [Part I](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation) will be analysed with a series of bioinformatic tools within the conda `RNA1` environment. These tools are used in a predetermined order to evaluate and improve the quality of paired-reads of each dataset before the alignment-based quantification of gene expression takes place.  

In more details, **preprocessing** consists of quality control (QC) of raw reads and, depending on the results, the reads are trimmed and filtered by length and quality (Phred score) to retain only high quality reads with minimal adapter contamination. After these QC steps, the **secondary analysis** starts in which the quality-improved datasets are mapped (aligned) to the reference genome, followed by flagging of PCR duplicates and quantifying aligned pair-end reads.

The way to assign all these processes in a predetermined order, ensuring that each step in these processes is consistent across data and platforms is by drafting a bioinformatic pipeline. Such pipelines utilize a programming or workflow language, allowing these processes to be portable, (if possible) parallelizable, consistent, and interoperable. These pipelines can be implemented using a variety of workflow systems, for example:

- **Bash**: simple, transparent, ideal for small workflows. More details:  
    - [How to write a bash script](https://www.youtube.com/watch?v=F-gskSl4pwQ)
    
- **Nextflow** & **Snakemake**: reproducible, scalable, cloud‑ready, container‑friendly. More details:  
    - [An introduction to Nextflow](https://www.commonwl.org)  
    - [An introduction to Snakemake](https://www.youtube.com/watch?v=tUTcfoMQl98&t=136s)
- **CWL** (Common Workflow Language) — standardized, portable workflows across platforms. More details:  
    - <https://www.commonwl.org>

- **WDL** (Workflow Description Language) + **Cromwell** — used by Broad Institute; strong support for large genomics pipelines. More details:   
    - <https://github.com/broadinstitute/cromwell>  
    - <https://www.youtube.com/watch?v=w0IUd-x_9NU>

- **Galaxy**: GUI*‑based workflow system for non‑programmers who would like to learn bioinformatics
(*GUI: Graphical User Interface). More details:  
    - <https://www.youtube.com/watch?v=k6fTVIR4GME>

In this tutorial, we will implement the pipeline using **Bash** and **Nextflow**. **Bash** is ideal for learning the underlying commands and logic of each step. **Nextflow** adds reproducibility, scalability, and the ability to resume failed jobs—valuable skills for real-world research.  

Both pipelines will cover **preprocessing** and **secondary analysis** of the datasets.  

> [!IMPORTANT]  
> **By the end of Part II, you will have:**
> - Cleaned, trimmed FASTQ files ready for alignment
> - Aligned reads in BAM format, sorted and indexed
> - Duplicate-marked BAM files for accurate quantification
> - A raw count matrix (`raw_counts.txt`) ready for differential expression analysis
> - Experience running the same pipeline with **Bash** and **Nextflow**



## Pipeline overview

The following table summarizes the steps, tools, inputs, and outputs, and description of the bulk RNA-seq pipeline implemented in this tutorial:


| **Step** | **Tool** | **Input** | **Output** | **Description** |
| :--- | :--- | :--- | :--- | :--- |
| 1. QC (Raw) | `FastQC` + `MultiQC` | Raw FASTQ files | QC reports (HTML + ZIP) | Assess raw read quality, GC content, adapter contamination, and overrepresented sequences |
| 2. Trimming | `Cutadapt` | Raw FASTQ files | Trimmed FASTQ (`.fastq.gz`) | Remove adapter sequences, trim low-quality bases, and filter reads by length |
| 3. QC (Trimmed) | `FastQC` + `MultiQC` | Trimmed FASTQ files | QC reports (HTML + ZIP) | Re-evaluate read quality after trimming to confirm improvement |
| 4. Alignment | `HISAT2` + `samtools sort` | Trimmed FASTQ files | Sorted BAM (`.sorted.bam`) | Align trimmed reads to the reference genome (GRCh38) and sort the resulting BAM files |
| 5. Duplicate Marking | `Picard MarkDuplicates` | Sorted BAM | Dedup BAM (`.dedup.bam`) + metrics | Flag PCR duplicates in aligned BAM files (without removing them, as required for RNA-seq) |
| 5.5. BAM Indexing | `samtools index` | Dedup BAM (`.dedup.bam`) | BAM index file (`.dedup.bam.bai`) | Create an index for the deduplicated BAM file to enable fast random access for downstream tools and visualization |
| 6. Strandedness | `RSeQC (infer_experiment.py)` | Dedup BAM + BED12 | Strandedness report (`.txt`) | Determine library strandedness to set the correct `-s` parameter for quantification |
| 7. Quantification | `featureCounts` | Dedup BAM + GTF | Raw count matrix (`raw_counts.txt`) | Count paired-end reads mapping to genes to generate a raw count matrix |
| 8. Post-Alignment QC | `RSeQC` + `MultiQC` | Dedup BAM + BED12 | QC reports + MultiQC summary | Assess alignment quality, read distribution, and splice junction annotation |
| 9. Visualization | `IGV` | Dedup BAM + BAI | Interactive genome browser view | Visualize aligned reads, splice junctions, and coverage across genomic regions |


Before the start of **prepro

(and take a decision whether to trim away bad quality bases and adapter sequences, as well as filter out reads based on bp length)




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
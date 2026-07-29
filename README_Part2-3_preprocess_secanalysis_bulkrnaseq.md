# Part II – Preprocessing & Secondary analysis

## Table of Contents

- [Introduction](#introduction)  
- [Pipeline overview](#pipeline-overview)  
- [Bash: Preprocessing](#bash-preprocessing)  
- [Bash: Alignment and mark duplicates](#bash-alignment-and-mark-duplicates)  
- [Bash: Gene-level paired-end read quantification](#iii-bash-gene-level-paired-end-read-quantification)  
- [Nextflow: Preprocessing](#i-nextflow-preprocessing)  
- [Nextflow: Alignment and mark duplicates](#ii-nextflow-alignment-and-mark-duplicates)  
- [Nextflow: Gene-level paired-end read quantification](#iii-nextflow-gene-level-paired-end-read-quantification)  
- [Visualization](#visualization)


## Introduction

In this second part of bulk RNA-seq analysis, the two datasets downloaded in [Part I](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation) will be analysed with a series of bioinformatic tools within the conda `RNA1` environment. These tools are used in a predetermined order to evaluate and improve the quality of paired-reads of each dataset before the alignment-based quantification of gene expression takes place.  

In more details, **preprocessing** consists of quality control (QC) of raw reads and, depending on the results, the reads are trimmed and filtered by length and quality (Phred score) to retain only high quality reads with minimal adapter contamination. After these QC steps, the **secondary analysis** starts with the mapping or alignment of the quality-improved datasets to the reference genome, followed by flagging of PCR duplicates and quantification of aligned pair-end reads.

The way to assign all these processes in a predetermined order, ensuring that each step in these processes is consistent across data and platforms is by drafting a ***bioinformatic pipeline***. Such pipelines utilize a programming or workflow language, allowing these processes to be portable, (if possible) parallelizable, consistent, and interoperable. These pipelines can be implemented using a variety of workflow systems, for example:

- **Bash**: simple, transparent, ideal for small workflows. More details:  
    - [How to write a bash script](https://www.youtube.com/watch?v=F-gskSl4pwQ)
    
- **Nextflow** & **Snakemake**: reproducible, scalable, cloud‑ready, container‑friendly. More details:  
    - [An introduction to Nextflow](https://www.youtube.com/watch?v=Gq1KiMJyNB4&t=1705s)  
    - [An introduction to Snakemake](https://www.youtube.com/watch?v=tUTcfoMQl98&t=136s)
- **CWL** (Common Workflow Language) — standardized, portable workflows across platforms. More details:  
    - <https://www.commonwl.org>

- **WDL** (Workflow Description Language) + **Cromwell** — used by Broad Institute; strong support for large genomics pipelines. More details:   
    - [Welcome to Cromwel - GitHub](https://github.com/broadinstitute/cromwell)  
    - [Intro to learning miniwdl for WDL](https://www.youtube.com/watch?v=w0IUd-x_9NU)

- **Galaxy**: GUI*‑based workflow system for non‑programmers who would like to learn bioinformatics
(*GUI: Graphical User Interface). More details:  
    - [What is Galaxy?](https://www.youtube.com/watch?v=k6fTVIR4GME)

In this tutorial, we will implement the pipeline using **Bash** and **Nextflow**. **Bash** is ideal for learning the underlying commands, bioinformatics tools, and the logic of each step. **Nextflow** adds reproducibility, scalability, and the ability to resume failed jobs, which are valuable skills for real-world research.  

Both pipelines will cover **preprocessing** and **secondary analysis** of the datasets.  

> [!IMPORTANT]  
> **By the end of Part II, you will have:**
> - Cleaned, trimmed FASTQ files ready for alignment
> - Aligned reads in BAM format, sorted and indexed
> - Duplicate-marked BAM files for accurate quantification
> - A raw count matrix (`raw_counts.txt`) ready for differential expression analysis
> - Experience running the same pipeline with **Bash** and **Nextflow**

---

## Pipeline overview

The following table summarizes the steps, tools, inputs, and outputs, and description of the bulk RNA-seq pipeline implemented in this tutorial:


| **Step** | **Tool** | **Input** | **Output** | **Description** |
| :--- | :--- | :--- | :--- | :--- |
| 1. QC (Raw) | `FastQC` + `MultiQC` | Raw FASTQ files | QC reports (HTML + ZIP) | Assess raw read quality, GC content, adapter contamination, and overrepresented sequences |
| 2. Trimming | `Cutadapt` | Raw FASTQ files | Trimmed FASTQ (`.fastq.gz`) | Remove adapter sequences, trim low-quality bases, and filter reads by length |
| 3. QC (Trimmed) | `FastQC` + `MultiQC` | Trimmed FASTQ files | QC reports (HTML + ZIP) | Re-evaluate read quality after trimming to confirm improvement |
| 4. Alignment | `HISAT2` +<br>`samtools sort` | Trimmed FASTQ files | Sorted BAM (`.sorted.bam`) | Align trimmed reads to the reference genome (GRCh38) and sort the resulting BAM files by genomic coordinates (chr and position) |
| 5. Duplicate Marking | `Picard MarkDuplicates` | Sorted BAM | Dedup BAM (`.dedup.bam`) + metrics | Flag PCR duplicates in aligned BAM files (without removing them, as required for RNA-seq) |
| 5.5. BAM Indexing | `samtools index` | Dedup BAM (`.dedup.bam`) | BAM index file (`.dedup.bam.bai`) | Create an index for the deduplicated BAM file to enable fast random access for downstream tools and visualization |
| 6. Strandedness | `RSeQC (infer_experiment.py)` | Dedup BAM + BED12 | Strandedness report (`.txt`) | Determine library strandedness to set the correct `-s` parameter for accurate **gene-level** quantification |
| 7. Quantification | `featureCounts` | Dedup BAM + GTF | Raw count matrix (`raw_counts.txt`) | Count paired-end reads mapping to genes to generate a **gene-level** raw count matrix for DE analysis|
| 8. Post-Alignment QC | `RSeQC` + `MultiQC` | Dedup BAM + BED12 | QC reports + MultiQC summary | Assess alignment quality, read distribution, and splice junction annotation |
| 9. Visualization | `IGV` | Dedup BAM + BAI | Interactive genome browser view | Visualize aligned reads, splice junctions, and coverage across genomic regions |

---

## Bash: Preprocessing  

### 1. QC report of raw datasets: FastQC and MultiQC  

As shown in [Part I](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation), to develop a bash script, you have to create an executable `.sh` file and next copy/paste/save the bash script below into the `.sh` file. Then, run it from `~/Bulk_rnaseq/scripts`.
<br>
**Steps**:  
1. Navigate to `Bulk_rnaseq/scripts`  
2. Create `RNA1_01_bulkrnaseq_preprocessing.sh`  
3. Grant execute permissions  

```bash
# Run these commands one by one
cd path/to/Bulk_rnaseq/scripts
touch RNA1_01_bulkrnaseq_preprocessing.sh
chmod u+x RNA1_01_bulkrnaseq_preprocessing.sh
```

4. Open the `.sh`. Use a text/script editor, e.g. nano, vim, etc. and copy/paste/save the **bash script** below   

**Bash script: raw datasets QC**  

```bash
#!/bin/bash

set -euo pipefail

# Set variables as path
DATA_DIR="$1"               # /path/to/Bulk_rnaseq
PROJECT="PRJNA437330"
PROJECT_PATH="$DATA_DIR/data/$PROJECT"
THREADS=4
RESULTS="$DATA_DIR/results"
QC_DIR="$RESULTS/qc_raw"
QC_DIR_FASTQC="$QC_DIR/fastq_raw"
QC_DIR_MULTIQC="$QC_DIR/multiqc_raw"
QC_DIR_FASTQC_TRIM="$RESULTS/qc_trimmed/fastq_trimmed"
QC_DIR_MULTIQC_TRIM="$RESULTS/qc_trimmed/multiqc_trimmed"
RAW_FASTQ_DIR=$PROJECT_PATH/*/raw_fastq     # To expand the '*' (placeholder for "SRR..." datasets) do not use quotation marks
TRIMMED="$RESULTS/trimmed"
LOGS="$RESULTS/logs"


# ------- QC fastq files -------

# Create a QC folder for raw fastq files
mkdir -p "$QC_DIR_FASTQC"
mkdir -p "$QC_DIR_MULTIQC"

echo "####################"
echo "## Running FASTQC ##"
echo "####################"

for fastq in $RAW_FASTQ_DIR/*.fastq.gz; do
    fastqc \
        --threads "$THREADS" \
        --outdir "$QC_DIR_FASTQC" \
        "$fastq"
done

# MultiQC report from raw fastq files
  echo "#####################"
  echo "## Running MultiQC ##"
  echo "#####################"

multiqc \
  "$QC_DIR_FASTQC" \
  -o "$QC_DIR_MULTIQC"    


# ------- Trimming & filtering -------            # to be continue
```

5. Running the script: In `Bulk_rnaseq/scripts` run the `.sh` with the working directory `~/Bulk_rnaseq`

```bash
./RNA1_01_bulkrnaseq_preprocessing.sh /path/to/Bulk_rnaseq
```

<br>

> [!NOTE]  
> The line `DATA_DIR="$1"` at the top of the bash script indicates the script where the path to the project folder directory, in this case is `/path/to/Bulk_rnaseq`, is located. The `$1` is a **command‑line argument** that you must provide when running the script.  
>
> **This is important**: You cannot simply run `./RNA1_01_bulkrnaseq_preprocessing.sh` because the script expects you to **pass the path** as an argument.
> 
> **Correct usage**:
> ```bash
> ./RNA1_01_bulkrnaseq_preprocessing.sh /path/to/Bulk_rnaseq
> ```
>
> **Incorrect usage** (will fail):
> ```bash
> ./RNA1_01_bulkrnaseq_preprocessing.sh
> ```
>
> Replace `/path/to/Bulk_rnaseq` with the actual path to your project folder.

<br>

6. Folder structure: Output files from **raw datasets QC**. See `~Bulk_rnaseq/results/qc_raw`

```bash
Bulk_rnaseq/
├── data
│   ├── PRJNA437330
│   │   ├── SRR6815993
│   │   │   └── raw_fastq
│   │   │       ├── SRR6815993_1.fastq.gz
│   │   │       └── SRR6815993_2.fastq.gz
│   │   └── SRR6816017
│   │       └── raw_fastq
│   │           ├── SRR6816017_1.fastq.gz
│   │           └── SRR6816017_2.fastq.gz
│   └── sra_PRJNA437330.sh
├── reference
├── results
│   ├── qc_raw
│   │   ├── fastq_raw
│   │   │   ├── SRR6815993_{1,2}_fastqc.html
│   │   │   ├── SRR6815993_{1,2}_fastqc.zip
│   │   │   ├── SRR6816017_{1,2}_fastqc.html
│   │   │   ├── SRR6816017_{1,2}_fastqc.zip
│   │   └── multiqc_raw
│   │       ├── multiqc_data
│   │       └── multiqc_report.html
└── scripts
    └── RNA1_01_bulkrnaseq_preprocessing.sh
```

7. **FastQC** and **MultiQC** reports

**FastQC** reports, for both samples and for each R1 and R2, show:

- Sequence length of 75bp and zero sequences flagged as poor quality  
- Good **Per base sequence quality**, but the first 5bp show lower quality than the rest  
- Warning sign in **Per base sequence content**  
- Levels of **duplication** high
- There are **overrepresented sequences**

<br>

**MultiQC** report: since the quality of reads and bp is excellent, the most important issue are the overrepresented sequences and the adapter content 

![**MultiQC of raw datasets**](images/multiqc_raw_samples_1.png)   



### 2. Trimming/filtering of reads + QC: Cutdapt + FastQC & MultiQC

For the second part of the **preprocessing** bash script:  

1. Copy/paste/save the **trimming/filtering of reads** and **post trimming QC** part to the `RNA1_01_bulkrnaseq_preprocessing.sh` file  

**Bash script: Trimming/filtering + QC**  
  
```bash
# ------- Trimming & filtering -------

# Create a trimming folder
mkdir -p "$TRIMMED"
mkdir -p "$LOGS"

# For looping each sample (SRR accession), process both R1 and R2
for SAMPLE_DIR in $PROJECT_PATH/*; do
  # Extract the sample name from the directory path (e.g., SRR6815993). 'basename' strips directory path and returns only the last component
  SAMPLE=$(basename "$SAMPLE_DIR")

  echo "######################"
  echo "## Running Cutadapt ##"
  echo "## Sample: $SAMPLE  ##"
  echo "######################"

  echo "$SAMPLE_DIR"

  # Define input and output file paths
  R1_IN="$SAMPLE_DIR/raw_fastq/${SAMPLE}_1.fastq.gz"
  R2_IN="$SAMPLE_DIR/raw_fastq/${SAMPLE}_2.fastq.gz"
  R1_OUT="$TRIMMED/${SAMPLE}_R1.trimmed.fastq.gz"
  R2_OUT="$TRIMMED/${SAMPLE}_R2.trimmed.fastq.gz"

  cutadapt \
    -j "$THREADS" \
    -u 5 -U 5 \
    -q 24,24 \
    -m 30 \
    --poly-a \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -b GTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
    -b TATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
    -b GGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
    -o "$R1_OUT" \
    -p "$R2_OUT" \
    "$R1_IN" "$R2_IN" \
    > "$LOGS/cutadapt_${SAMPLE}.log" 2>&1

    echo "✅ Trimming complete for $SAMPLE"
    echo "   R1 output: $R1_OUT"
    echo "   R2 output: $R2_OUT"
    echo

done

# ------- QC fastq files -------

# Create a QC folder for raw fastq files
mkdir -p "$QC_DIR_FASTQC_TRIM"
mkdir -p "$QC_DIR_MULTIQC_TRIM"

# No looping for FASTQC this time because the trimmed files are all in $TRIMMED
echo "####################"
echo "## Running FASTQC ##"
echo "####################"

fastqc \
  --threads "$THREADS" \
  --outdir "$QC_DIR_FASTQC_TRIM" \
  $TRIMMED/*.trimmed.fastq.gz


# MultiQC report from raw fastq files
echo "#####################"
echo "## Running MultiQC ##"
echo "#####################"

multiqc \
  $QC_DIR_FASTQC_TRIM \
  -o $QC_DIR_MULTIQC_TRIM

```
<br>
 
 
> [!NOTE]  
> **How the for loop works**:
> 
> The `for` loop processes both samples (`SRR6815993` and `SRR6816017`) automatically. 
>
> There is no single correct way to write a `for` loop in Bash. You can implement it in different ways, as long as it correctly iterates over the input files.
>
> `SAMPLE_DIR` contains the full paths within `$PROJECT_PATH/*`, which are:
>
> - `path/to/Bulk_rnaseq/data/PRJNA437330/SRR6815993`  
> - `path/to/Bulk_rnaseq/data/PRJNA437330/SRR6816017`  
>
> `SAMPLE` contains only the **sample ID**, extracted from `SAMPLE_DIR` using the `basename` command:  
> - `SRR6815993`  
> - `SRR6816017`  
>
> **Cutadapt options used**:  
>
> `-a / -A`: Remove Illumina Nextera adapter sequences from Read 1 and Read 2.  
> `--poly-a`: Trim poly-A tails from Read 1 and poly-T heads from Read 2 after adapter trimming.  
> `-u / -U`: Remove the first 5 bases from Read 1 and Read 2 before quality trimming.  
> `-q 24,24`: Trim low-quality bases from both the 5′ and 3′ ends of Read 1 and Read 2 using a Phred quality cutoff of 24. Quality trimming is performed before adapter trimming.  
> `-m`: Discard reads shorter than the specified minimum length after trimming.  
> `"$LOGS/cutadapt_${SAMPLE}.log" 2>&1`: Redirect both the standard output (stdout) and standard error (stderr) to a single log file for each sample. During a successful run, the log mainly contains Cutadapt's processing summary; if any errors occur, they are written to the same file.  
>
> **Order of Cutadapt operations**: 
> 1. Remove the first 5 bases (`-u / -U`).
> 2. Perform quality trimming (`-q`).
> 3. Remove adapter sequences (`-a / -A`).
> 4. Trim poly-A tails (R1) and poly-T heads (R2) (`--poly-a`).
> 5. Discard reads shorter than 30 nt (`-m`).





2. Folder structure: Output files from trimmed datasets and post QC.  

See:
- `~/Bulk_rnaseq/results/qc_trimmed`  
- `~/Bulk_rnaseq/results/trimmed`  
- `~/Bulk_rnaseq/results/logs`

```bash
Bulk_rnaseq/
├── data
│   ├── PRJNA437330
│   │   ├── SRR6815993
│   │   │   └── raw_fastq
│   │   │       ├── SRR6815993_1.fastq.gz
│   │   │       └── SRR6815993_2.fastq.gz
│   │   └── SRR6816017
│   │       └── raw_fastq
│   │           ├── SRR6816017_1.fastq.gz
│   │           └── SRR6816017_2.fastq.gz
│   └── sra_PRJNA437330.sh
├── reference
├── results
│   ├── logs
│   │   ├── cutadapt_SRR6815993.log
│   │   └── cutadapt_SRR6816017.log
│   ├── qc_raw
│   ├── qc_trimmed
│   │   ├── fastq_trimmed
│   │   │   ├── SRR6815993_R1.trimmed_fastqc.{html,zip}
│   │   │   ├── SRR6815993_R2.trimmed_fastqc.{html,zip}
│   │   │   ├── SRR6816017_R1.trimmed_fastqc.{html,zip}
│   │   │   └── SRR6816017_R2.trimmed_fastqc.{html,zip}
│   │   └── multiqc_trimmed
│   │       ├── multiqc_data
│   │       └── multiqc_report.html
│   └── trimmed
│       ├── SRR6815993_R1.trimmed.fastq.gz
│       ├── SRR6815993_R2.trimmed.fastq.gz
│       ├── SRR6816017_R1.trimmed.fastq.gz
│       └── SRR6816017_R2.trimmed.fastq.gz
└── scripts
    └── RNA1_01_bulkrnaseq_preprocessing.sh

```

3. **FastQC** and **MultiQC** reports

**FastQC** reports, for both samples and for each R1 and R2, show:

- Sequence length of 75bp and zero sequences flagged as poor quality  
- Good **Per base sequence quality**, but the first 5bp show lower quality than the rest  
- Warning sign in **Per base sequence content**  
- Levels of **duplication** high
- There are **overrepresented sequences**

<br>

**MultiQC** report: since the quality of reads and bp is excellent, the most important issue are the overrepresented sequences and the adapter content 

![**MultiQC of raw datasets**](images/multiqc_raw_samples_2.png)   

<br>


<br>
<br>
If you have reached the end of **PART I**, I congratulate you!!  
Continue to the 👉 [Part II – Secondary analysis](README_Part2-3_secondary_bulkrnaseq.md), where you'll start with the preprocessing analysis to alignment till the generation of raw counts tables, using bash and nextflow scripting explained step-by-step.

Go back to the top of 👉 [Part I – Setup & data preparation](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation)

Go to the main page 👉 [Bulk RNA-seq Tutorial](README.md)
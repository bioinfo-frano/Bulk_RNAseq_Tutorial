# Part II – Preprocessing & Secondary analysis

## Table of Contents

- [Introduction](#introduction)  
- [Pipeline overview](#pipeline-overview)  
- [Bash: Preprocessing](#bash-preprocessing)  
    - [1. QC report of raw datasets](#bash-preprocessing)  
    - [2. Trimming / filtering of reads & QC](#2-trimmingfiltering--qc-cutadapt--fastqc--multiqc)
- [Bash: Secondary analysis](#bash-secondary-analysis) 
    - [1. Alignment / marking duplicates & QC](#bash-secondary-analysis)  
    - [2. Determination of strandedness](#2-determination-of-strandedness)  
    - [featureCounts: Gene-level paired-end read quantification](#3-featurecounts-gene-level-paired-end-read-quantification)  
- [Nextflow: Preprocessing](#i-nextflow-preprocessing)  
- [Nextflow: Alignment and mark duplicates](#ii-nextflow-alignment-and-mark-duplicates)  
- [Nextflow: Gene-level paired-end read quantification](#iii-nextflow-gene-level-paired-end-read-quantification)  
- [Visualization](#visualization)


## Introduction

In this second part of bulk RNA-seq analysis, the two datasets downloaded in [Part I](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation) will be analysed with a series of bioinformatic tools within the conda `RNA1` environment. These tools are used in a predetermined order to evaluate and improve the quality of paired-reads of each dataset before the alignment-based quantification of gene expression takes place.  

In more details, the **preprocessing** consists of quality control (QC) of raw reads and, depending on the results, the reads are trimmed and filtered by length and quality, according to their Phred score, to retain only high quality reads with minimal adapter contamination. After these QC steps, the **secondary analysis** starts with the mapping or alignment of the quality-improved datasets to the reference genome, followed by flagging of PCR duplicates and quantification of aligned pair-end reads.

The way to assign all these processes in a predetermined order, ensuring that each step in these processes is consistent across data and platforms is by drafting a ***bioinformatic pipeline***. Such pipelines utilize a programming or workflow language and bioinformatic tools, allowing these processes to be portable, (if possible) parallelizable, consistent, and interoperable. These pipelines can be implemented using a variety of workflow systems, for example:

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

- **Galaxy**: *GUI‑based workflow system for non‑programmers who would like to learn bioinformatics. More details:  
    - [What is Galaxy?](https://www.youtube.com/watch?v=k6fTVIR4GME)  
    - ***GUI**: Graphical User Interface  

In this tutorial, we will implement a pipeline using **Bash** and **Nextflow**. **Bash** is ideal for learning the underlying commands, bioinformatics tools, and the logic of each step. **Nextflow** adds reproducibility, scalability, and the ability to resume failed jobs, which are valuable skills for real-world research.  

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

### 1. QC report of raw datasets: FastQC & MultiQC  

As shown in [Part I - Find & download paired-end RNA-seq datasets](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation), to develop a bash script, you have to create an executable `.sh` file and next copy/paste/save the bash script below into the `.sh` file. Then, run it from `~/Bulk_rnaseq/scripts`.
<br>
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

6. **Folder structure**: Output files from **raw datasets QC**. See `~Bulk_rnaseq/results/qc_raw`

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
- There are **overrepresented sequences** and observed **adapter content**

<br>

**MultiQC** report:  
Since the quality of reads and bp is excellent, the most important issue are the overrepresented sequences and the adapter content. See the snapshot of the MutiQC report below, showing the overrepresented and adapter sequences content in R1 & R2 of both samples (`SRR6816017`, `SRR6815993`).

<br>

![**MultiQC of raw datasets**](images/multiqc_raw_samples_1.png)   



### 2. Trimming/filtering & QC: Cutadapt & FastQC & MultiQC

For the second part of the **preprocessing** bash script.
<br>
<br>
**Steps**:  
1. Copy/paste/save the **trimming & filtering of reads** and **post trimming QC** part to the `RNA1_01_bulkrnaseq_preprocessing.sh` file 

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



2. **Folder structure**: Output files from trimmed datasets and post QC.  

See:
- `~/Bulk_rnaseq/results/qc_trimmed`  
- `~/Bulk_rnaseq/results/trimmed`  
- `~/Bulk_rnaseq/results/logs`
- `~/Bulk_rnaseq/scripts`

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

**FastQC** reports after trimming and filtering of reads, for both samples and for each R1 and R2, show:

- Sequence length of 30 - 70bp and zero sequences flagged as poor quality  
- Improved **Per base sequence quality**
- Warning sign in **Per base sequence content**  
- Levels of **duplication** high
- Less than 1% of reads with **overrepresented sequences** and less than 0.1% with **adapter content**

<br>

**MultiQC** report:  

- **Sequence Length Distribution**: There's a wider range of sequence lengths; however, most of the reads are 70bp
- **Mean Quality Scores**: Quality of reads improved even more
- **Overrepresented sequences by sample** & **Adapter Content**: content of overrepresented sequences and adapters is almost negligible. See the snapshot of the MutiQC report below, showing the overrepresented and adapter sequences content in R1 & R2 of both samples (`SRR6816017`, `SRR6815993`) after trimming.

<br>

![**MultiQC of raw datasets**](images/multiqc_aftertrimming_samples_1.png)   

---


## Bash: Secondary analysis

### 1. Alignment, mark duplicates and QC: HISAT2 + MarkDuplicates + MultiQC

1. Create another `.sh` script to perform read alignment with **HISAT2**, mark potential PCR duplicates, and generate a post-alignment QC report.

```bash
# Run these commands one by one
cd path/to/Bulk_rnaseq/scripts
touch RNA1_02_bulkrnaseq_alignment_markdup.sh
chmod u+x RNA1_02_bulkrnaseq_alignment_markdup.sh
```

2. Open the `.sh`. Use a text/script editor, e.g. nano, vim, etc. and copy/paste/save the **bash script** below   

**Bash script: Alignment + mark duplicates + QC**  

```bash
#!/bin/bash

set -euo pipefail

# Set variables as path
DATA_DIR="$1"               # /path/to/Bulk_rnaseq
PROJECT="PRJNA437330"
PROJECT_PATH="$DATA_DIR/data/$PROJECT"
THREADS=4
RESULTS="$DATA_DIR/results"
RAW_FASTQ_DIR=$PROJECT_PATH/*/raw_fastq     # To expand the '*' (placeholder for "SRR..." datasets) do not use quotation marks
TRIMMED="$RESULTS/trimmed"
LOGS="$RESULTS/logs"
HISAT2_INDEX="$DATA_DIR/reference/hisat2_index/grch38_tran"
ALIGNMENT="$RESULTS/alignment"
QC_POST_ALIGN="$RESULTS/qc_post_align"

# ------- Aligment: HISAT2 -------

# Create folders for alignment (in case they don't exist)
mkdir -p "$ALIGNMENT"
mkdir -p "$LOGS"

# Define sample IDs and names as indexed arrays (compatible with Bash 3.x)
SAMPLES=("SRR6815993" "SRR6816017")
SAMPLE_NAMES=("6h_Mock" "6h_STM-D23580_inv")

# For looping the alignment per sample
for i in "${!SAMPLES[@]}"; do
  SAMPLE_ID="${SAMPLES[$i]}"
  SAMPLE_NAME="${SAMPLE_NAMES[$i]}"

  # Input trimmed fastq files
  R1_TRIM="$TRIMMED/${SAMPLE_ID}_R1.trimmed.fastq.gz"
  R2_TRIM="$TRIMMED/${SAMPLE_ID}_R2.trimmed.fastq.gz"

  echo "########################"
  echo "## Running HISAT2     ##"
  echo "## Sample: $SAMPLE_ID ##"
  echo "########################"

  hisat2 \
    -x "$HISAT2_INDEX/genome_tran" \
    -1 "$R1_TRIM" \
    -2 "$R2_TRIM" \
    --rg-id "${SAMPLE_ID}" \
    --rg "SM:${SAMPLE_NAME}" \
    --rg "LB:${SAMPLE_ID}" \
    --rg "PL:ILLUMINA" \
    --rg "PU:unknown" \
    --new-summary \
    --summary-file "$LOGS/${SAMPLE_ID}.hisat2.log" \
    -p "$THREADS" \
  | samtools sort \
    -@ "$THREADS" \
    -m 500M \
    -o "$ALIGNMENT/${SAMPLE_ID}.sorted.bam"

    echo "✅ Alignment complete for $SAMPLE_ID"
    echo "Output: $ALIGNMENT/${SAMPLE_ID}.sorted.bam"
    echo

done

# ------- Mark duplicates: Picard -------

# For looping the Markduplicates per sample
for i in "${!SAMPLES[@]}"; do
  SAMPLE_ID="${SAMPLES[$i]}"

  echo "################################"
  echo "## Running MarkDuplicates     ##"
  echo "## Sample: $SAMPLE_ID         ##"
  echo "################################"

  picard MarkDuplicates \
    I="$ALIGNMENT/${SAMPLE_ID}.sorted.bam" \
    O="$ALIGNMENT/${SAMPLE_ID}.dedup.bam" \
    M="$LOGS/${SAMPLE_ID}_dedup_metrics.txt" \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=SILENT

    echo "✅ Duplicate marking complete for $SAMPLE_ID"
    echo "Output: $ALIGNMENT/${SAMPLE_ID}.dedup.bam"
    echo "Metrics: $LOGS/${SAMPLE_ID}_dedup_metrics.txt"

    # Index the dedup BAM for IGV visualization
    echo "🔍 Indexing BAM file..."
    samtools index "$ALIGNMENT/${SAMPLE_ID}.dedup.bam"
    echo "✅ Indexing complete"

done

# ------- MultiQC Post-alignment -------

mkdir -p "$QC_POST_ALIGN"

echo "################################"
echo "## Running MultiQC            ##"
echo "## Directory: $QC_POST_ALIGN  ##"
echo "################################"

multiqc \
  "$LOGS" \
  -o "$QC_POST_ALIGN"

```
    
<br>
 
 
> [!NOTE]  
> **How the alignment and mark duplicates loops work**:
>
> At this stage of the pipeline, the input files are the trimmed paired-end FASTQ files located in `~/Bulk_rnaseq/results/trimmed`   
>
> The `for` loop processes both samples (`SRR6815993` and `SRR6816017`) automatically. The script defines a list called `SAMPLES`, which stores the sequencing accession IDs. A second list, `SAMPLE_NAMES`, stores descriptive sample names that are added to the BAM file as **read group** (**RG**) metadata during alignment. During each iteration, the loop aligns one sample, adds the corresponding **RG** information, sorts the alignments with `samtools sort`, and saves a separate HISAT2 log file for that sample.
>
> Including **RG** information during alignment (by aligners such as **HISAT2**, **BWA-MEM**, and **STAR**), is considered good practice because it allows downstream tools to distinguish sequencing libraries and samples. The **RG** provides information about:
>
> - **ID**: Sample identifier (`SRR6815993` and `SRR6816017`)
> - **SM**: Name of biological sample (`6h_Mock`, `6h_STM-D23580_inv`)
> - **LB**: Library. The authors stated that "Barcoded Illumina sequencing libraries (Nextera XT...) were generated...", which means that samples `SRR6815993` and `SRR6816017` had unique barcode identifiers (each barcoded sample represents a separate physical library). Therefore, `LB` is set to the sample ID `SRR6815993` and `SRR6816017`
> - **PL**: Platform information (`ILLUMINA`)
> - **PU**: Platform unit. A platform unit should identify the flowcell + lane + index, e.g. `HF7K2DMXX.1.ATCACG`. This information is not available in the SRA metadata and is stripped from the FASTQ headers. Thus, `PU` is set to the value `unknown`. Picard MarkDuplicates does not require this field (`PU`).
>
> A second `for` loop runs **Picard** `MarkDuplicates` on each sorted BAM file and flags those duplicated reads, generating a `.dedup.bam` file per sample, and one duplication metrics report `*_dedup_metrics.txt` per sample.  
>
> Duplicate reads are flagged, **NOT removed**, because of the option: `REMOVE_DUPLICATES=false` is used. This preserves all reads while allowing downstream tools to identify PCR duplicates if needed.
>
> Finally, each `.dedup.bam` file is indexed with `samtools index`, generating a corresponding `.dedup.bam.bai` file. The `.dedup.bam` and `.dedup.bam.bai` are required for efficient read alignment visualization with **IGV**.  

    
3. **Folder structure**: Output files alignment and post QC.  

See:
- `~/Bulk_rnaseq/results/alignment`  
- `~/Bulk_rnaseq/results/logs`  
- `~/Bulk_rnaseq/results/qc_post_align`
- `~/Bulk_rnaseq/scripts`

```bash
Bulk_rnaseq/
├── data
├── reference
├── results
│   ├── alignment
│   │   ├── SRR6815993.dedup.bam
│   │   ├── SRR6815993.dedup.bam.bai
│   │   ├── SRR6815993.sorted.bam
│   │   ├── SRR6815993_chr_prefix.txt
│   │   ├── SRR6816017.dedup.bam
│   │   ├── SRR6816017.dedup.bam.bai
│   │   ├── SRR6816017.sorted.bam
│   │   └── SRR6816017_chr_prefix.txt
│   ├── logs
│   │   ├── SRR6815993.hisat2.log
│   │   ├── SRR6815993_dedup_metrics.txt
│   │   ├── SRR6816017.hisat2.log
│   │   ├── SRR6816017_dedup_metrics.txt
│   │   ├── cutadapt_SRR6815993.log
│   │   └── cutadapt_SRR6816017.log
│   ├── qc_post_align
│   │   ├── multiqc_data
│   │   └── multiqc_report.html
│   ├── qc_raw
│   ├── qc_trimmed
│   └── trimmed
└── scripts
    ├── RNA1_01_bulkrnaseq_preprocessing.sh
    └── RNA1_02_bulkrnaseq_alignment_markdup.sh
```
    
4. **MultiQC** report  

- **HISAT2**: Pair-ends (PE) reads mapped uniquely  
  - `SRR6815993`: 83.1% ✅  
  - `SRR6816017`: 77.4% ✅  
  
- **Mark Duplicates**:   

| **Sample**    | **Unique Pairs** | **Duplicate Pairs nonoptical** |
| :---          | :---             | :---                           | 
| `SRR6815993`  | 64.6%            | 26.6%                          |
| `SRR6816017`  | 53.7%            | 37.1%                          |

<br>

> [!IMPORTANT]  
> It's importat to notice the difference between **optical** and **non optical duplicates** and 
> **Optical duplicates**: It's a type of sequencing artefact, which happens when a single amplification cluster is incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument 👉 [Documentation: MarkDuplicates (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360036834611-MarkDuplicates-Picard).   
> To calculate this parameter, the optical duplicates parameter was not calculated. To do so, it's necessary to add the option `--READ_NAME_REGEX` and `--OPTICAL_DUPLICATE_PIXEL_DISTANCE ` to **MarkDuplicates** chunk.  
>
> **Duplicate Pairs nonoptical**: Another type of artifact, where identical DNA/RNA fragments are generated during library preparation, primarily via PCR amplification, rather than from optical sensor errors on the sequencer.
>
> <https://www.reddit.com/r/bioinformatics/comments/1as62s3/can_anyone_help_me_figure_out_what_is_duplicate/>

<br>

![**MultiQC of HISAT2 and MarkDuplicates**](images/multiqc_hisat2_md_samples_1.png)   

<br>

> [!NOTE]  
> **Optical Duplicates in SRA Data**:
> 
> In the MarkDuplicates metrics file, you may see `READ_PAIR_OPTICAL_DUPLICATES = 0`. This is because SRA datasets often have **stripped read names**, which means that flow cell metadata doesn't have information about tile, cluster, and X/Y coordinates. Without this information, Picard cannot distinguish optical duplicates (flow cell artifacts) from PCR duplicates.
> 
> **How to check your data**:
> ```bash
> zcat SRR6815993_1.fastq.gz | head -1
> ```
> 
> **If you see a header like this one**:
> ```
> @SRR6815993.1 1 length=75
> ```
>
> This means that the header doesn't have the read names, it lacks flow cell/tile/cluster info → **Optical duplicates cannot be calculated**.
> 
> **If you see**:
> ```
> @A00489:123:HF7K2DMXX:1:1101:10000:10000 1:N:0:ATCACG
> ```
> The read names contain full flow cell metadata → Optical duplicates **can** be calculated.
> 
> **What to do**:
> - If your bulk RA-seq data is from SRA (like this tutorial): **Accept that optical duplicates cannot be calculable** – this is normal.
> - If your data is from your own sequencing run: Use the original FASTQ files with full read names.
> - Even without optical duplicate detection, other duplication metrics (`PERCENT_DUPLICATION`, `ESTIMATED_LIBRARY_SIZE`, and the duplicate set histogram) are still valid and useful for QC.
> 
> **Key takeaway**: The absence of optical duplicate detection is **not an error** – it's a limitation of the SRA data format. It does not affect the quality of your gene-level count matrix.

<br>

- **Cutadapt**: 
  - Filtered Reads
    - `SRR6815993`: 97.2% ✅  
    - `SRR6816017`: 96.2% ✅  
  - Trimmed Sequence Lengths (3'): shows the amount of reads trimmed in x amount of bases from their 3' end. Example, when x-axis shows 10bp length trimmmed and 5000 read counts in y-axis means that there are 5000 reads had 10bp trimmed from 3' end.

<br>


![**MultiQC of HISAT2 and MarkDuplicates**](images/multiqc_cutadapt_samples_1.png)  


<br>

### 2. Determination of strandedness

During transcription, RNA polymerase reads the **template** strand in the 3'→5' direction and synthesizes mRNA in the 5'→3' direction. The opposite strand is the **coding** strand, which has the same sequence as the mRNA. Genes can be encoded on either strand, producing overlapping mRNAs when transcribed from opposite strands.
During the library preparation, using a **stranded** (strand-specific) library retains the information about the original strand of DNA from which the mRNA was transcribed. This improves **transcript annotation** and **quantification**, particularly when distinguishing overlapping genes, antisense transcripts, and non-coding RNAs transcribed from opposite strands. In contrast, conventional **non-stranded** RNA-seq libraries lose information about the strand of origin during double-stranded cDNA preparation.

See these papers for more details:  
- [Comprehensive comparative analysis of strand-specific RNA sequencing methods](https://www.nature.com/articles/nmeth.1491)  
- [Comparison of stranded and non-stranded RNA-seq transcriptome profiling and investigation of gene overlap](https://link.springer.com/article/10.1186/s12864-015-1876-7)  
- [Signal & Kahlke, 2021: how_are_we_stranded_here: quick determination of RNA‑Seq strandedness](https://pmc.ncbi.nlm.nih.gov/articles/PMC8783475/)  

**Assuming a stranded library as unstranded** can result in **over 10% false positives** and **over 6% false negatives** in downstream differential expression results (Signal et al., *BMC Bioinformatics*, 2022).  
The strandedness information **is not available** for RNA-sequencing samples in repositories such as ENA or SRA, and **publications often do not report this information in the methods**. In fact, a randomised investigation of 50 ENA paired-end studies found that only 56% explicitly stated or mentioned strandedness in their methods (Signal et al., *BMC Bioinformatics*, 2022). Therefore, **it is important to determine the strandedness of our datasets**.

`infer_experiment.py` from the package **RSeQC** is one of the tools used to determine the strandedness of RNA-seq data. The tool requires a **.bed** and a **.bam** alignment file. It compares the orientation of aligned reads against known gene annotations to infer strandedness.
  
`infer_experiment.py` reports three possible outcomes:

- **Stranded** (forward): reads map to the same strand as the transcript

- **Reversely stranded** (reverse): reads map to the opposite strand

- **Unstranded**: reads map to both strands with roughly equal frequency


**Before** testing strandedness, you must verify that the chromosome naming between your HISAT2 alignment file, **.bam**, and the annotation file, **BED12**, match.  

If you don't remember how the **BED12** file was created, check the **chapter V** "Create a BED12 file" from 👉 [Part I – Setup & data preparation](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation)

<br> Let's observe both files

**SRR6815993.dedup.bam**

```bash
samtools view -h SRR6815993.dedup.bam | head -26
```

Output:

```bash
@HD	VN:1.6	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:11	LN:135086622
...
@SQ	SN:20	LN:64444167
@SQ	SN:21	LN:46709983
@SQ	SN:3	LN:198295559
...
@SQ	SN:8	LN:145138636
@SQ	SN:9	LN:138394717
@SQ	SN:MT	LN:16569
@SQ	SN:X	LN:156040895
@SQ	SN:Y	LN:57227415
```

**SRR6816017.dedup.bam**

```bash
samtools view -h SRR6816017.dedup.bam | head -26
```

Output:

```bash
@HD	VN:1.6	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:11	LN:135086622
...
@SQ	SN:20	LN:64444167
@SQ	SN:21	LN:46709983
@SQ	SN:3	LN:198295559
...
@SQ	SN:8	LN:145138636
@SQ	SN:9	LN:138394717
@SQ	SN:MT	LN:16569
@SQ	SN:X	LN:156040895
@SQ	SN:Y	LN:57227415
```

Chromosomes don't show prefix `chr` for both datasets.

Since it was created a **BED12** file without `chr` prefix called `gencode.v38.annotation.nochr.clean.bed`, then we can use this for the strandedness analysis.

<br>

Now, check the strandedness:

```bash
STRANDED="$RESULTS/strandedness"
BED12_NOCHR="$DATA_DIR/reference/intervals/gencode.v38.annotation.nochr.clean.bed"
COUNTS_DIR="$RESULTS/raw_counts"

# ------- Strandedness -------

mkdir -p "$STRANDED"

for i in "${!SAMPLES[@]}"; do
  SAMPLE_ID="${SAMPLES[$i]}"

  echo "###########################################"
  echo "## Running RSeQC (infer_experiment.py)   ##"
  echo "## Sample: $SAMPLE_ID                    ##"
  echo "###########################################"

infer_experiment.py \
  -r "$BED12_NOCHR" \
  -i "$ALIGNMENT/${SAMPLE_ID}.dedup.bam" \
  > "$STRANDED/${SAMPLE_ID}_strandedness.txt"

done
```

<br>

Output:  
Sample: **SRR6815993.dedup.bam**  


```bash
This is PairEnd Data
Fraction of reads failed to determine: 0.1107
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4421
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4472
```
<br>

Output:  
Sample: **SRR6816017.dedup.bam**  


```bash
This is PairEnd Data
Fraction of reads failed to determine: 0.1150
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4428
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4422
```
<br>

**Interpretation**  

Read:  
- R1 = 1  
- R2 = 2  
Read strand: + or -  
Gene strand: + or -  

So `1++` means: "Read 1 mapped to the + strand, and the gene is on the + strand."  
**Positive / Sense strand**: forward or coding strand that shares the same sequence direction and 5' - 3' orientation as the corresponding mRNA.

Group/Pattern 1: `"1++,1--,2+-,2-+": 0.4428` → **Forward stranded**   
Group/Pattern 2: `"1+-, 1-+, 2++, 2--": 0.4422` → **Reverse stranded**  

Both configurations occur at almost exactly the same frequency: 44.28% vs 44.22%

This is typical of an **unstranded** library

<br>

| **Library type**     | **Fraction 1** (`1++,1--,2+-,2-+`) | **Fraction 2** (`1+-,1-+,2++,2--`) |
|:---------------------|:-----------------------------------|:-----------------------------------|
| **Forward stranded** | **High** (> 0.7)                   | Low (< 0.2)                        |
| **Reverse stranded** | Low (< 0.2)                        | **High** (> 0.7)                   |
| **Unstranded**       | **~0.45**                          | **~0.45**                          |


<br>

The strandedness **describes** the relationship between the sequencing read and the original mRNA.

| Library type | Relationship | What `infer_experiment.py` shows |
| :--- | :--- | :--- |
| Forward stranded (sense) | R1 maps to the same strand as the mRNA | Group 1 high (`1++,1--,2+-,2-+`) |
| Reverse stranded (antisense) | R1 maps to the opposite strand of the mRNA | Group 2 high (`1+-,1-+,2++,2--`) |
| Unstranded | No strand information preserved | Both groups ≈ 45% |

<br>

> [!IMPORTANT]  
> A stranded kit preserves strand information for every transcript, **regardless of whether** the gene sits on the + strand or the − strand of the chromosome. The kit does not "prefer" genes on one strand. Whether the strandedness of reads is forward, reverse, or unstranded is crucially important for setting up the right `-s` parameter in featureCounts. If your library is reverse-stranded but you use `-s 1` (which means **forward stranded**), many reads will be discarded — sometimes half of them. This is the most common cause **of missing or under-counted genes.**

A kit is either unstranded or stranded, and if it's stranded, it is specifically either forward-stranded or reverse-stranded. The forward/reverse distinction is not optional — it's baked into the chemistry of the kit.

The distinction comes from which strand of the cDNA is sequenced. Two common mechanisms:
| Mechanism | Result | Example kits |
| :--- | :--- | :--- |
| **dUTP method**: dUTP is incorporated during second-strand synthesis, then the second strand is degraded before sequencing. Only the first strand (antisense to mRNA) is read. | Reverse stranded | Illumina TruSeq Stranded mRNA, Illumina Stranded Total RNA, NEB Ultra II Directional |
| **Ligation-based directional methods**: Adapters are ligated in a specific orientation that preserves the original mRNA strand as the read. | Forward stranded | Some older ligation-based kits, certain Lexogen protocols |

The **dUTP method is by far the most common library kit in modern RNA-seq**, which is why most modern stranded kits are reverse-stranded.

The kit's chemistry determines the read orientation regardless of which gene you're looking at. So:

- A **reverse-stranded kit** always produces reads where R1 is antisense to the mRNA, for every gene, whether it's on the + or − strand. In other words, in a reverse-stranded library, R1 is the **reverse complement of the mRNA (antisense)**
- A **forward-stranded kit** always produces reads where R1 is sense to the mRNA, for every gene, whether it's on the + or − strand. In other words, in a forward-stranded library, R1 is the **same sequence as the mRNA (sense)**

This is why `infer_experiment.py` reports a single strandedness value for the whole library, not per-gene. The strandedness is a property of the protocol, not of any individual gene.

<br>

**Practical example**

**Reverse-Stranded Example: Two Genes, Same Kit**  
Take a **reverse-stranded kit** (dUTP method — the most common modern protocol) and two genes:  

| Gene | Gene strand | R1 maps to | R2 maps to | R1 vs. mRNA | `infer_experiment.py` codes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Gene A | + | − | + | Antisense | R1: `1-+` · R2: `2++` |
| Gene B | − | + | − | Antisense | R1: `1+-` · R2: `2--` |

All four codes fall into **Group 2** (`1+-,1-+,2++,2--`), which is diagnostic of a **reverse-stranded** library.  

In `infer_experiment.py` output, this group would show a fraction **> 0.7**, while Group 1 (`1++,1--,2+-,2-+`) would be **< 0.2**:  

```bash
This is PairEnd Data
Fraction of reads failed to determine: 0.11
Fraction of reads explained by "1++,1--,2+-,2-+": 0.10   ← Low (Group 1)
Fraction of reads explained by "1+-,1-+,2++,2--": 0.79   ← High (Group 2) → reverse stranded
```

**Forward-Stranded Example: Two Genes, Same Kit**  
Take a **forward-stranded kit** (ligation-based directional method) and the same two genes:  

| Gene | Gene strand | R1 maps to | R2 maps to | R1 vs. mRNA | `infer_experiment.py` codes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Gene A | + | + | - | Sense | R1: `1++` · R2: `2+-` |
| Gene B | − | - | + | Sense | R1: `1--` · R2: `2-+` |


All four codes fall into **Group 1** (`1++,1--,2+-,2-+`), which is diagnostic of a **forward-stranded** library.  

In `infer_experiment.py` output, this group would show a fraction > 0.7, while Group 2 (1+-,1-+,2++,2--) would be < 0.2:

```bash
This is PairEnd Data
Fraction of reads failed to determine: 0.11
Fraction of reads explained by "1++,1--,2+-,2-+": 0.79   ← High (Group 1) → forward stranded
Fraction of reads explained by "1+-,1-+,2++,2--": 0.10   ← Low (Group 2)
```

**Why is this important?**: (because) The strandedness information must be passed to `featureCounts` via the `-s` parameter to ensure correct gene counting.

| Library type | Group 1 (`1++,1--,2+-,2-+`) | Group 2 (`1+-,1-+,2++,2--`) | featureCounts `-s` parameter |
| :--- | :--- | :--- | :--- |
| Forward stranded | High (> 0.7) | Low (< 0.2) | `-s 1` |
| Reverse stranded | Low (< 0.2) | High (> 0.7) | `-s 2` |
| Unstranded | ~0.45 | ~0.45 | `-s 0` |

<br>

**In summary**: "Stranded" refers to the protocol's ability to preserve strand-of-origin information, not to a preference for one chromosomal strand over the other. Both forward-stranded and reverse-stranded kits work for all genes. The difference is purely in the read orientation relative to the mRNA — and that difference is what `infer_experiment.py` detects and what `featureCounts -s` must match. Within stranded kits, the read orientation can be forward (sense) or reverse (antisense), depending on the chemistry. Most modern Illumina stranded kits are reverse-stranded (dUTP method). The strandedness is a property of the protocol, not of any individual gene.

<br>

### 3. featureCounts: Gene-level paired-end read quantification

Up to this point, the reads have been trimmed and quality checked, aligned, and the strandedness determined. The next step is the counting of "genes". In this case, a '**gene**' refers to a genomic locus, and the '**count**' represents the number of read pairs (fragments) that align to its exons. Each counted fragment (R1 + R2) originates from an mRNA (or non‑coding RNA) molecule present in the sample, so its abundance serves as a proxy for gene expression.  
The tool that calculates this, which considers the strandedness and the metadata of genes, is **featureCounts**.  
<br>
**featureCounts** is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads (Subread website, see documentation below).

**Documentation**  
1. [Subread](https://subread.sourceforge.net)  
2. [featureCounts](https://subread.sourceforge.net/featureCounts.html)  

<br>


**Bash script: featureCounts**

```bash
BED12_NOCHR="$DATA_DIR/reference/intervals/gencode.v38.annotation.nochr.clean.bed"
COUNTS_DIR="$RESULTS/raw_counts"
INTERVAL_GTF="$DATA_DIR/reference/intervals/gencode.v38.annotation.gtf.gz"
QC_POST_STRAND_COUNTS="$RESULTS/qc_strandedness_rawcounts"
QC_POST_ALIGN_RSEQC="$RESULTS/qc_post_align_rseqc"

mkdir -p "$COUNTS_DIR"

# Build an array of BAM files
BAM_FILES=()   # ← Initialize empty array

for i in "${!SAMPLES[@]}"; do
  SAMPLE_ID="${SAMPLES[$i]}"
  BAM_FILES+=("$ALIGNMENT/${SAMPLE_ID}.dedup.bam")
done

echo "BAM files to process: ${BAM_FILES[@]}"
echo ""

echo "####################################"
echo "## Running featureCounts          ##"
echo "## Generating combined count      ##"
echo "## matrix for all samples         ##"
echo "####################################"

echo

featureCounts \
  -T "$THREADS" \
  --countReadPairs \
  # -s 0: unstranded library (see strandedness section)
  -s 0 \
  -a "$INTERVAL_GTF" \
  -o "$COUNTS_DIR/raw_counts.txt" \
  -t exon \
  -g gene_id \
  --extraAttributes gene_name \
  -p \
  "${BAM_FILES[@]}"

echo
echo "########################################"
echo "## featureCounts completed successfully"
echo "########################################"
echo "Output: $COUNTS_DIR/raw_counts.txt"
echo "Summary: $COUNTS_DIR/raw_counts.txt.summary"

echo
echo "Summary statistics:"
echo "========================================"
cat "$COUNTS_DIR/raw_counts.txt.summary"
echo "========================================"

```

<br>

The output of featureCounts will show a `.txt` file containing a table, showing metadata per gene in columns 1-6, and the so-called **raw counts** on the seventh column, which shows the counting information. **The count column is named after each BAM file** and the columns are: Geneid, Chr (chromosome), Start, End, Strand, Length, SRR6815993.dedup.bam (dataset) 



<br>
<br>
<br>

```bash
**MAKE A GTF FILE WITHOUT "chr" PREFIX, AND THE RUN THE FEATURE COUNT AND QC ALL OVER AGAIN**
# 1. Confirm the GTF now uses "MT" for mitochondria (not "M" or "chrM")
zcat gencode.v38.annotation.nochr.gtf.gz | grep -v '^#' | awk '$1=="MT"' | head -1

# 2. Confirm the output Chr column shows "1" not "chr1"
head -3 results/raw_counts/raw_counts.txt | cut -f1,2

# Create a GTF without "chr" prefix, matching the BAM
zcat gencode.v38.annotation.gtf.gz | sed 's/^chr//' | gzip > gencode.v38.annotation.nochr.gtf.gz
INTERVAL_GTF="$DATA_DIR/reference/intervals/gencode.v38.annotation.nochr.gtf.gz"
```

<br>
<br>

If you have reached the end of **PART I**, I congratulate you!!  
Continue to the 👉 [Part II – Secondary analysis](README_Part2-3_secondary_bulkrnaseq.md), where you'll start with the preprocessing analysis to alignment till the generation of raw counts tables, using bash and nextflow scripting explained step-by-step.

Go back to the top of 👉 [Part I – Setup & data preparation](README_Part1-3_setup_bulkrnaseq.md#part-i--setup--data-preparation)

Go to the main page 👉 [Bulk RNA-seq Tutorial](README.md)
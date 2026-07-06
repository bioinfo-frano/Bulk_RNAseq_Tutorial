# Part I – Setup & data preparation

## Table of Contents

- [Introduction](#introduction)
- [Creating the computing environment for bulk RNAseq analysis](#creating-the-computing-environment-for-bulk-rnaseq-analysis)
    - [I. Create a folder structure](#i-create-a-folder-structure)
    - [II. Create a Conda environment](#ii-create-a-conda-environment)
    - [III. Find & download FASTQ datasets](#iii-find--download-fastq-datasets)
    - [IV. Download a pre-built HISAT2 genome index](#iv-download-a-pre-built-hisat2-genome-index)



## Introduction

The analysis of transcriptomic datasets independently can be quite difficult if you don't have the proper guidance and, importantly, enough patience, time and computational resources. At the end, you would have to send your datasets to an external bioinformatician or try other options all of which imply financial costs. This is the logical solution when the statistics, tables and plots are urgently needed for the submission of scientific manuscripts or when preparing seminars.  
The purpose of this tutorial is to show that you can independently analyse your data relying on a personal computer (e.g., laptop) or workstation, which has limited computational resources.  
For the sake of learning, how to analyse your RNAseq datasets, ideally, you should have some basic knowledge on command line, bash scripting, R programming, and python. However, if you don't have it, don't worry, **learn by doing it!**

I would strongly suggest the following tutorials, so that you train yourself in these topics.
  
- [Learn Mac Terminal Basics - macmostvideo - YouTube](https://www.youtube.com/watch?v=ZkoEHvG3GI8)   
- [Command Line Basics for Beginners - Full Course - freeCodeCamp.org - YouTube](https://www.youtube.com/watch?v=mABpAI-pCw0)    
- [MASTERING Command Prompt Basics! | Tutorial (for Windows) - Skill Foundry - YouTube](https://www.youtube.com/watch?v=QBWX_4ho8D4)   
- [Bash Scripting Tutorial for Beginners  - TechWorld with Nana - YouTube](https://www.youtube.com/watch?v=PNhq_4d-5ek)   

Optional, but highly recommended too:   
  
- [R programming in one hour - a crash course for beginners - R Programming 101 - YouTube](https://www.youtube.com/watch?v=eR-XRSKsuR4&t=176s)   
- [Python Full Course for Beginners - Programming with Mosh - YouTube](https://www.youtube.com/watch?v=K5KVEU3aaeQ)   
  
Before starting, I would also recommend you to read the [Part I - Preparation & setup -  Introduction](https://github.com/bioinfo-frano/NGS_Workflow_Tutorial/blob/main/README_Part1-3_setup.md), where you can read a bit more about cloud computing alternatives and on what **FASTQ** files are.

With that foundational knowledge in mind, let's now set up our local environment for the actual analysis.  
  
From all analyses in this tutorial, the **alignment step** is the most <u>**computational demanding**</u>. Therefore, since we are limited in terms of computational power, this tutorial will provide pipelines for the analysis of small numbers of RNA datasets that can be processed comfortably on standard workstations or laptops. Later on, when working in **R**, it will be possible to expand the amount of RNA datasets by downloading a pre-aligned raw counts. Let's start.


## Creating the computing environment for bulk RNAseq analysis

I.	Create a folder structure  

II.	Create conda environment

III.	Find & download FASTQ datasets from a published scientific paper

IV. Download a pre-built HISAT2 genome indexes (e.g., Homo sapiens GRCh38/hg38)

> [!NOTE]   
> This guide was developed and tested on macOS running on Intel processors. Users on Apple Silicon (M1/M2/…/M5) or Linux systems may need to adapt certain steps.


> [!IMPORTANT] **Conda prerequisites:**  
> This guide assumes that Miniconda3 is already installed on your computer. If not, please consult the official [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/macos.html) or watch this [YouTube](https://www.youtube.com/watch?v=OH0E7FIHyQo) video.

When Miniconda is already installed, you should see the `(base)` environment activated in your Terminal.

## I. Create a folder structure  


## I. Create a specific conda environment called `DNA`

### a) Create the environment with these dependencies.

> **IMPORTANT**
You don't need to deactivate from `(base)` when creating a new Conda environment.

```bash
conda create -n DNA \
  -c conda-forge -c bioconda -c defaults \
  python=3.9 \
  openjdk=17 \
  perl=5.32 \
  fastqc \
  multiqc \
  cutadapt \
  bwa \
  samtools \
  picard \
  htslib \
  gatk4 \
  bcftools \
  vcftools \
  snpeff \
  ensembl-vep \
  bedtools \
  coreutils \
  pigz \
  pbzip2 \
  pandas \
  numpy \
  matplotlib \
  seaborn \
  -y
```

### b) Sanity checks - post installation

  1. List Conda environments (the new `DNA` environment should appear):

  ```bash
  conda env list
  ```

  2. Activate the new environment.

```bash
  conda activate DNA
```

  3. Verify dependencies in `DNA` and version:

```bash
  conda list
```

  gatk --version
  
  java -version
  
  bwa
  
  samtools --version
  
  vep --help
  
  snpEff -version
  

You may find the following error when running `snpEff -version`: 

```bash
Error: LinkageError
UnsupportedClassVersionError
class file version 65.0
Java runtime only recognizes up to 61.0
```
**Explanation**: 
snpEff 5.3 was compiled with Java 21
You are running Java 17 (`openjdk=17`)
Java 17 cannot run Java 21 bytecode

**Recommendation (for stability)**: 
Downgrade snpEff:

```bash
conda install -n DNA -c bioconda snpeff=5.1
```

Expected output:

```bash
SnpEff	5.1d	2022-04-19
```

### c) Reproducibility of conda `DNA`

To export the environment, run in the Terminal:

```bash
conda env export --no-builds > DNA_conda_environment_full.yml
```

This generates a `DNA_conda_environment_full.yml` file containing all software dependencies.

>**Note**: 
The `--no-builds` flag is important because it removes platform-specific build strings that can break installation on different systems.

To reproduce this environment on another system, run:

```bash
conda env create -f DNA_conda_environment_full.yml
conda activate DNA
conda list
```

Finally, `conda list` will display the installed packages and confirm that the `DNA` conda environment was created successfully.

The **.yml** file is available here 👉 [DNA_conda_environment_full.yml](DNA_conda_environment_full.yml)

---

## II. Create the folder structure

All FASTQ files, reference (e.g. human) genome and scripts should be located into specific folders. Below is a recommended folder structure:

```bash
Genomics_cancer/
├── reference/                 # Reference genomes and known sites
│   └── GRCh38/
│       ├── fasta/             # Reference FASTA files
│       └── known_sites/       # Known variant sites (e.g. dbSNP, Mills)
│       └── bed/               # Genomic interval files (.bed) 
│       └── somatic_resources/ # Population and somatic reference.
├── data/
│   └── SRA_ID/                # Sample-specific directory (e.g. SRR30536566)
│       ├── raw_fastq/         # Original FASTQ files
│       ├── qc/                # FastQC / MultiQC reports
│       ├── trimmed/           # Adapter- and quality-trimmed FASTQ files
│       ├── aligned/           # BAM files and indexes
│       ├── variants/          # VCF files
│       └── annotation/        # Annotated variants
├── scripts/                   # scripts
└── logs/                      # Log files from pipeline execution

bed/  
# Genomic interval files (.bed) defining target regions (e.g. gene panels,
# exons, amplicons). Used for read filtering, coverage calculation,
# and restricting variant calling to clinically relevant regions.

somatic_resources/  
# Population and somatic reference resources required for somatic variant calling,
# especially for GATK Mutect2 (e.g. Panel of Normals, gnomAD allele frequencies,
# germline resource VCFs).
```

Multiple samples can be processed by creating one directory per SRA accession under `data/`.

In Terminal, create all directories at once::

```bash
mkdir -p Genomics_cancer/{reference/GRCh38/{fasta,known_sites},data/SRR30536566/{raw_fastq,qc,trimmed,aligned,variants,annotation},scripts,logs}
```

---

## III. Find & download small-sized FASTQ datasets for cancer gene panels

Downloading FASTQ files directly from the SRA web interface is **not recommended**, because:

1. R1 and R2 reads may be merged into a single file

2. There is no guarantee that the FASTQ files represent original raw reads (they may be reconstructed from alignments)

Instead, use **SRA-Tools**.

### 1. Install SRA-Tools in a separate Conda environment

```bash
conda create -n sra \
>   -c conda-forge \
>   -c bioconda \
>   python=3.10 \
>   sra-tools=3
```

Activate `(sra)` environment:

```bash
conda activate sra
```

Verify installation

```bash
fasterq-dump --version
```

Expected output:

```bash
fasterq-dump : 3.2.1			# If output is 'fasterq-dump : 2.9.6', SRA-tools won't work. Update!
```

> [!NOTE]
> The next step is configuring the **SRA Toolkit** in order to "***access public and, optionally, controlled-access data in the cloud***"    
> To do so, follow the steps in <https://www.uvm.edu/vacc/docs/beyond_basics/sratoolkit/>, which tells you how to configure the **SRA Toolkit** via the interactive menu.    
> The interactive menu can be reached with the command `vdb-config --interactive` or `vdb-config -i` through Terminal. Here, navigate through **SRA configuration** to set up, for example, the cache directory (default is `~/ncbi/public/sra/`).  
>
> By configuring the cache, tools like `prefetch` will store the downloaded `*.sra` files in `~/ncbi/public/sra/`. Later, `fasterq-dump` can read from the cache to split into R1 and R2 FASTQ files.  
> The cache can be cleared with `cache-mgr -c` or `cache-mgr --clear`. Use with caution, as this removes all cached `.sra` files.   
>
> If this configuration is **not** set up, the `~/ncbi/public/sra/` won't exist or be used, and `cache-mgr` will be appear disabled. Despite this, `prefetch`, `vdb-validate`, `fasterq-dump` and other tools from **SRA Toolkit** will still work. However, each time you use `prefetch`, the `.sra` files will be placed in your **current working directory** instead of the cache (`~/ncbi/public/sra/`).  



### 2. Find small FASTQ files from [SRA](https://www.ncbi.nlm.nih.gov/sra)

> [!IMPORTANT]  
> **Controlled-access human genomic data (dbGaP):**  
> Many sequencing datasets derived from human subjects (e.g., cancer or germline studies) are considered sensitive because they may contain identifiable genetic information.  
>
> As a result:
> - Raw sequencing data (FASTQ, BAM, CRAM) are often **not publicly available**  
> - Access is controlled through the database of Genotypes and Phenotypes (**dbGaP**) managed by the **NCBI**  
>
> To obtain access, researchers typically must:
> - Be affiliated with a recognized institution
> - Be a senior investigator/scientist/clinician, tenure-track investigator or staff scientist.
> - Submit a formal data access request  
> - Agree to data use limitations  
> - Provide **Institutional Review Board (IRB)** approval for some datasets
>
> Therefore, fully open-access human tumor-only or tumor–normal datasets are limited.  
> For training purposes, use:
> - Public SRA/ENA datasets explicitly marked as open-access (**CONSENT PUBLIC**)
> - Small targeted panel (**Assay Type: Targeted-Capture**) datasets when available  

**Suggestion:**
- Use keywords: **targeted, illumina, cancer, genomic, Homo sapiens**
- Example dataset: `SRR30536566`

### 3. Inspect dataset structure (size, source, format) before downloading

```bash
vdb-dump --info SRR30536566
```
Expected output:

```bash
acc    : SRR30536566
path   : https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30536566/SRR30536566
remote : https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30536566/SRR30536566
size   : 339,709,019
type   : Database
platf  : SRA_PLATFORM_ILLUMINA
SEQ    : 3,892,036
SCHEMA : NCBI:SRA:Illumina:db#2
TIME   : 0x0000000066d7e0e7 (09/04/2024 06:24)
FMT    : sharq
FMTVER : 3.0.11
LDR    : general-loader.3.0.8
LDRVER : 3.0.8
LDRDATE: Sep 11 2023 (9/11/2023 0:0)
```
**Key fields to inspect:**
```bash
size   : 339,709,019                       # ~340 MB dataset. Compressed SRA size
SCHEMA : NCBI:SRA:Illumina:db#2            # Illumina reads forward and reverse
FMT    : sharq                             # Compressed format
```

**Table 3: Interpretation of** `vdb-dump --info`.

| **SCHEMA**                       | **FMT** | **What it really is**                              | **Suitable for FASTQ-first pipelines?** |
| -------------------------------- | ------- | -------------------------------------------------- | --------------------------------------- |
| `NCBI:SRA:GenericFastq`          | `FASTQ` | Raw FASTQ (platform-agnostic)                      | ✅ Yes                                   |
| `NCBI:SRA:GenericFastq`          | `sharq` | Raw FASTQ stored in SRA compressed format          | ✅ Yes                                   |
| `NCBI:SRA:Illumina`              | `FASTQ` | Raw FASTQ (Illumina-native representation)         | ✅ Yes                                   |
| `NCBI:SRA:Illumina`              | `sharq` | Raw Illumina FASTQ stored in compressed SRA format | ✅ Yes                                   |
| `NCBI:align:db:alignment_sorted` | `FASTQ` | FASTQ reconstructed from aligned reads             | ⚠️ Acceptable (with caveats)             |
| `NCBI:align:db:alignment_sorted` | `BAM`   | Internally stored aligned reads (BAM-like)         | ⚠️ Acceptable via FASTQ reconstruction   |
| `NCBI:align:db:alignment_sorted` | `CRAM`  | Internally stored aligned reads (CRAM-like)        | ⚠️ Acceptable via FASTQ reconstruction   |

If `SCHEMA` contains `align`, the dataset is **not raw**, and FASTQ files will be **reconstructed from aligned reads**.


**Table 4: Example of SRA dataset with raw and aligned FASTQ files.**

| Accession   | SCHEMA       | FMT   | Raw FASTQ? |
| ----------- | ------------ | ----- | ---------- |
| ERR12140864 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR20701732 | align        | FASTQ | ❌ No       | 
| SRR35865210 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR35529667 | SRA:Illumina | sharq | ✅ Yes      | 
| SRR32679397 | SRA:Illumina | sharq | ✅ Yes      | 

As shown in Table 3 and Table 4, both raw and aligned FASTQ datasets can be used but the `SRR20701732` is not really a raw FASTQ file. However, when downloading `SRR20701732` via `fasterq-dump`, this dataset will be converted to FASTQ.

> [!NOTE]
> Some SRA datasets (especially targeted panels or clinical studies) are stored internally as aligned reads (`SCHEMA: NCBI:align:db:alignment_sorted`).
>
> When downloading with `fasterq-dump` (or via **ENA**), reads are reconstructed into standard FASTQ format.
>
>These reconstructed FASTQ files:
>
> - preserve sequences, quality scores, and pairing
>
> - are suitable for QC, alignment, and somatic variant calling (including tumor–normal analysis)
>
> However:
>
> - they may differ slightly from the original sequencer output (e.g., read order or filtering)
>
> - therefore, truly raw datasets (`SCHEMA: NCBI:SRA`) are preferred when available

### 4. Download the SRA dataset

- Go to working directory `~/data/SRR30536566`

```bash
fasterq-dump SRR30536566 \    
  --split-files \		      	  
  --threads 4 \		        	  
  --outdir raw_fastq	  	  	
```
- `fasterq-dump SRR30536566` → SRA Run accession ID
- `--split-files`            → Splits paired-end reads into two separate FASTQ files: SRR30536566_1.fastq → forward reads | SRR30536566_2.fastq → reverse reads
- `--threads 4`              → Uses 4 CPU threads for parallel processing.
- `--outdir raw_fastq`       → Specifies the output directory where the FASTQ files will be saved. Alternatively, use the full absolute path and download the dataset from any directory.



Expected output:

```bash
spots read      : 3,892,036
reads read      : 7,784,072
reads written   : 7,784,072
```

> [!NOTE]  
> The `fasterq-dump` can be executed from any directory, but always use `--outdir raw_fastq` to maintain a clean and reproducible folder structure.

> [!IMPORTANT]
> **Do not get confuse**: SRA stores data in different formats (`FMT`) from raw FASTQ to aligned BAM/CRAM. However, even if `vdb-dump --info` reports `FMT: BAM` or `FMT: CRAM`, running `fasterq-dump` will still generate FASTQ files.
> Therefore:
> 
> - **You will not obtain BAM/CRAM files using** `fasterq-dump`
> - **You will always obtain FASTQ files suitable for downstream analysis**


Then, in folder `~/raw_fastq` there should be two fastq files, each having ~1.34 GB:

```bash
raw_fastq/
├── SRR30536566_1.fastq
└── SRR30536566_2.fastq
```

### 5. Compress FASTQ files

In `~/raw_fastq` directory:

```bash
gzip *.fastq
```

Expected output:

```bash
SRR30536566_1.fastq.gz
SRR30536566_2.fastq.gz
```

| Representation     | Approx size |
| ------------------ | ----------- |
| SRA (`sharq`)      | ~340 MB     |
| FASTQ (plain text) | ~2.6 GB (both) |
| FASTQ (gzipped)    | ~530 MB (both) |


Alternatively, compress during download with `fastq-dump` (but it's slow, not recommended for large data, probably deprecated)

```bash
fastq-dump SRR15506490 \
  --split-files \
  --threads 4 \
  --gzip \
  --outdir raw_fastq
```
Comparison between `fastq-dump` and `fasterq-dump`
| Tool           | Speed | Supports gzip directly | Recommendation         |
| -------------- | ----- | ---------------------- | ---------------------- |
| `fastq-dump`   | Slow  | ✅ Yes (`--gzip`)       | ❌ Avoid for large data |
| `fasterq-dump` | Fast  | ❌ No                   | ✅ Preferred tool       |



### 6. Verify FASTQ integrity

**Count reads**:

```bash
gzcat SRR30536566_1.fastq.gz | wc -l | awk '{print $1/4}'
```
Expected output:

`3892036`

```bash
gzcat SRR30536566_2.fastq.gz | wc -l | awk '{print $1/4}'
```
Expected output:

`3892036`

Inspect read structure:

```bash
zless SRR30536566_1.fastq.gz | head -n 8
```
**Expected output**:

```bash
@SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
TNTGTTTTTCCCTCCTGTTTTTTTTTTTTTTTCCTTAAAACCTACCATTTTTTCGGCATTTTGTTTTTTTTTTTTTTTTTTTTCCTTAGATTCATAATCA
+SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
A#AFFJJJJJJJJJJJJJJJJJJJJJJJJJJJ--<-<-77--7-7--<-7<7<77-7--<<----7A-7FF7F-<<7FAFF<--------7---7-777-
@SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
CNCTATTCTACCGGAAGCTGGAAGCAGCCTGACAGGATGTGTGGTGCCCAAGTCTGCACAGTGAGGTGGGGAGTGAGGGCCGCAGGCAGAGGGCAAGGGA
+SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
A#AFFFJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJJJJJAJJFJFJJJJJJJJJFJJJJJJJJJJJJJFJAFJJ7F<J
```

```bash
zless SRR30536566_1.fastq.gz| tail -n 8
```
Expected output:

```bash
@SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
GGAGGTGGGGCCCGGTGGAGGGTGATTGGATCATGGGGGTGGATTTCTCATCAGTGGTTTAGCACTATCCCCCTTAGTGCTGTGGTCACAATAGTGAGTG
+SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJJJ
@SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
ATGGTGATTGCATCTAATGTTTTCCTGTTATAGGGCAAATAATAGTGGTGATCTGGGTAATAGTTTCTCCAAATAATGACAAGCAGAAGTATACTCTGAA
+SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
```

```bash
zless SRR30536566_2.fastq.gz | head -n 8
```
Expected output:

```bash
@SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
NCACAAAAACAGCAGATGAAAGAGTCCTGGAAAGCTGCCTTCAAACTACCTCAGGGTTTTCACTAACTTTTACATAGACAGCTTTGATCTTGACCAGGAA
+SRR30536566.1 K00133:507:H2N2NBBXY:7:1101:1621:1068 length=100
#-AAFJ-F77J<<J<J-<FJJFJA--7-<FJJFJF-FA---7FJJF-FJA--F7A77AFJAJJ-AJJAA-7AAJ<J-FAF<7---7F-7----7--7)FF
@SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
NGCCTCTCCGCTGCCCTCTGCCCTGTGCCCCCTGCCCCCTGTCCCCTGTCCCTTGCCCTCTGCCTGCGGCCCTCACTCCCCACCTCACTGTGCAGACTTG
+SRR30536566.2 K00133:507:H2N2NBBXY:7:1101:2849:1086 length=100
#AAFFJFJJAJJAJJJJJJJJJJ-F<FJJJJJJ-AJJJJFJJJJJJF7FFJJJ7-FJJ<JFFJJ-FJFJJJJAJJJ7FJJJFFJAJFF7--<JF<AJAFF
```

```bash
zless SRR30536566_2.fastq.gz| tail -n 8
```
Expected output:

```bash
@SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
TGGTTCTTTGGCCTATATAGACTTCTGCTTTTGGTGGGGTGGGGGGCCAGGAAGCTTCCAATCATGGCAGAAGGCAAAAAGGGAGGAGTCATCTTACACA
+SRR30536566.3892035 K00133:507:H2N2NBBXY:7:2218:21339:47823 length=100
AAAFFJFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJ<7JJAAFFFJJJJJJFF77AFJF<JJJJJ
@SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
TATACTTGCCCTGATATTCTAAAACACAGAGTTTTAGTTGTTCAGAGGATAGCAACATACTTCGAGTTTTTTTCCTGATTGCTTCAGCAATTACTTGTTC
+SRR30536566.3892036 K00133:507:H2N2NBBXY:7:2218:23267:47823 length=100
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJ
```

Check:

- `@` header
- sequence line
- `+`
- quality scores

> [!NOTE]
> Using SRA Toolkit (`fasterq-dump`) may not always work. SRA datasets are often stored in a highly compressed internal format and may appear relatively small. However, converting them to FASTQ using `fasterq-dump` can be computationally intensive and requires **large temporary disk space**. As a result, `fasterq-dump` can be slow for big datasets or simple fail due to disk space limitations.

> [!TIP]
> When SRA Toolkit fails, download the FASTQ files directly from [ENA - European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home), specially for large datasets.
>
> **Protocol**:
>
> - Copy the dataset/run ID accession (**SRR#**) from SRA 
>
> - Paste it into the **ENA** "Search" bar
>
> - Locate the run in the results table
>
> - In the column "**Generated FASTQ files: FTP**" find and select both paired-end files
>
> - Click on either "**Get dowload script**", to retrieve a `wget` command to download the files or "**Download selected files**" directly.
>
> - Perform an **MD5 checksum** of each downloaded `*.fastq.gz` file and compare it with the value provided in ENA (column "**Generated FASTQ files: MD5**"). The checksums must match to ensure file integrity after download.
>
> Example: `md5 SRRXXXXXX_1.fastq.gz`
>
> Output: `MD5 (SRRXXXXXX_1.fastq.gz) = 21cfge23...`
>
> This is often a more efficient and reliable alternative than SRA, and it reduces disk usage.


> [!IMPORTANT]
> FASTQ files downloaded from **ENA** are often **reconstructed from aligned data** (if original was BAM) and, same as with `fasterq-dump`, they are perfectly valid for analysis.
> SRA is the authoritative source, but ENA is often the most practical way to obtain FASTQ files.

| Scenario                 | Recommended method |
| ------------------------ | ------------------ |
| Small datasets (<1 GB)   | SRA Toolkit        |
| Large datasets (>2–3 GB) | ENA + `wget`       |

> [!WARNING]
> Downloading FASTQ files from ENA is fast and convenient, but it is important to understand that these files **are not always identical to the original raw sequencing data stored in SRA**. In practice, ENA FASTQ files may differ from those obtained via `fasterq-dump` in several ways:
>
> - **Amount of read counts**: Some reads might be missing in **ENA** files.
> 
> - **Quality (Phread) scoremay be altered**: ENA files can contain binned or recalibrated Phred scores, leading to lower apparent base quality in tools like FastQC compared to SRA-derived FASTQ.
>
> **How to detect inconsistencies**
>
> To verify dataset integrity, always compare:
>
> - **Read counts (spots)**. Compare the number of reads reported in SRA with the downloaded dataset from **ENA**:
>
> `zcat sample.fastq.gz | wc -l | awk '{print $1/4}'`
>
> - **Quality profiles (FastQC)**. Check whether the per-base quality matches expectations from SRA metadata.
>
> - **Read length distribution**. A high proportion of short reads (<50 bp) in ENA data may indicate preprocessing or filtering.
>
> **Best practice** 
>
> Although ENA is faster and more convenient, the recommended standard for reproducible and high-integrity analysis is: 
>
> `fasterq-dump → gzip → downstream analysis`
> 
> This ensures that:
>
> - You are working with the original, unmodified sequencing data
>
> - Results are fully reproducible
> 
> - No hidden preprocessing steps affect your analysis


### 7. Using a Bash script for downloading more than one FASTQ file using SRA Toolkit, including evaluation of data integrity

```text
#!/bin/bash

set -euo pipefail

DATASETS=("SRR6815...1" "SRR6816...2" "SRR6816...3")  # EXAMPLES!!

for DATASET in "${DATASETS[@]}"; do
  echo
  echo "processing dataset: $DATASET"
  echo
  prefetch "$DATASET"
  vdb-validate "$DATASET" || { echo "Validation failed for $DATASET"; exit 1; }
  echo "downloading dataset: $DATASET"

  fasterq-dump "$DATASET" \
  --split-files \
  --threads 4 \
  --outdir "$DATASET"/raw_fastq

  echo "Dataset $DATASET downloaded successfully in $PWD"
  echo "Compressing"
  pigz -p 4 "$DATASET"/raw_fastq/*.fastq
  echo "Compression of $DATASET done!"
  echo "Removing $DATASET.sra file"
  rm -f "$DATASET/$DATASET.sra"  # If SRA Tools wasn't configured, then removed *.sra with 'rm -f', otherwise comment this out
  echo "Evaluating integrity"
  zgrep -c "@" "$DATASET"/raw_fastq/*.fastq.gz   # Alternative: zcat "$DATASET"/raw_fastq/*_1.fastq.gz | wc -l | awk '{print $1/4}'
  echo "Integrity evaluation finished"
  echo

done
```

| Command           | Function |  Benefit  |
| ----------------- | -------- |-------- |
| `prefetch`        | Downloads `.sra` file to cache        | Resume support, separates download from conversion   |
| `vdb-validate`    | Checks integrity of `.sra` file       | Prevents conversion of corrupted data   |
| `fasterq-dump`    | Converts `.sra` → FASTQ <br>(and can also directly download)   |Final step, generates usable FASTQ files |


---

## IV. Download a reference human genome (GRCh38) and indexes

The reference human genome is actually a bundle of files and it can be found on the GATK [website](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle), specifically in the link provided in the Resource Bundle hosted on [Google Cloud Buckets - gcp-public-data--broad-references ](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0)

In Google Cloud: ***Buckets/gcp-public-data--broad-references/hg38/v0*** is possible to find the whole reference human genome bundle (called **hg38** (informal name) or **GRCh38** (Genome Reference Consortium human build 38)). There you can find the following files:

- Homo_sapiens_assembly38.fasta 

- Homo_sapiens_assembly38.dict 

- Homo_sapiens_assembly38.fasta.fai

- Homo_sapiens_assembly38.fasta.64.alt

- Homo_sapiens_assembly38.fasta.64.amb

- Homo_sapiens_assembly38.fasta.64.ann

- Homo_sapiens_assembly38.fasta.64.bwt

- Homo_sapiens_assembly38.fasta.64.pac 

- Homo_sapiens_assembly38.fasta.64.sa

These are BWA index components, but now in 64-bit addressing mode.

| File      | Purpose                                   |
| --------- | ----------------------------------------- |
| `.64.amb` | Ambiguous bases (Ns)                      |
| `.64.ann` | Sequence names & lengths                  |
| `.64.bwt` | Burrows–Wheeler Transform                 |
| `.64.pac` | Packed reference sequence                 |
| `.64.sa`  | Suffix array (64-bit)                     |
| `.64.alt` | ALT contig ↔ primary contig relationships |


In Terminal, go to `/Genomics_cancer/reference/GRCh38/fasta`

**Download reference genome and indexes**: Run 👉 [0_wget_Hsapiens_assem38.sh](bash_scripts/0_wget_Hsapiens_assem38.sh)


# Reference Genome

**Genome:** GRCh38 / hg38  
**Source:** Broad Institute – GATK Resource Bundle  
**Bucket:** gcp-public-data--broad-references/hg38/v0

> **Why verify reference genome integrity?**
The reference genome FASTA and its associated index files (.fai and .dict) are foundational to all downstream NGS analyses. Any corruption or mismatch between these files can lead to alignment failures, GATK errors, or subtle coordinate inconsistencies that compromise variant calling results. Verifying file integrity using MD5 checksums ensures that the reference genome was downloaded correctly and that all associated index files correspond exactly to the same FASTA sequence, thereby guaranteeing reproducibility and analytical correctness.

MD5 checksums provide a fingerprint of a file.
If the checksum matches the expected value, the file is bit-by-bit identical to the original.

**Step 1: Go to the reference FASTA directory**

```bash
cd Genomics_cancer/reference/GRCh38/fasta
```

**Step 2: Compute MD5 checksums**

```bash
md5 Homo_sapiens_assembly38.fasta
md5 Homo_sapiens_assembly38.fasta.fai
md5 Homo_sapiens_assembly38.dict
```

Expected output:

```bash
7ff134953dcca8c8997453bbb80b6b5e  Homo_sapiens_assembly38.fasta
f76371b113734a56cde236bc0372de0a  Homo_sapiens_assembly38.fasta.fai
3884c62eb0e53fa92459ed9bff133ae6  Homo_sapiens_assembly38.dict
```


### (Final) Folder structure: before starting the analysis.

```bash
Genomics_cancer/
├── reference/                 
│   └── GRCh38/
│       ├── fasta/
│       │   ├── Homo_sapiens_assembly38.fasta
│       │   ├── Homo_sapiens_assembly38.fasta.fai
│       │   └── Homo_sapiens_assembly38.dict
│       │   └── Homo_sapiens_assembly38.fasta.64.amb     
│       │   └── Homo_sapiens_assembly38.fasta.64.ann     
│       │   └── Homo_sapiens_assembly38.fasta.64.bwt     
│       │   └── Homo_sapiens_assembly38.fasta.64.pac    
│       │   └── Homo_sapiens_assembly38.fasta.64.sa     
│       │   └── Homo_sapiens_assembly38.fasta.64.alt                  
│       └── known_sites/       
│       └── bed/               
│       └── somatic_resources/ 
├── data/
│   └── SRR30536566/                
│       ├── raw_fastq/
│       │   ├── SRR30536566_1.fastq.gz
│       │   └── SRR30536566_2.fastq.gz
│       ├── qc/                
│       ├── trimmed/           
│       ├── aligned/           
│       ├── variants/          
│       └── annotation/        
├── scripts/                   
└── logs/                      
```

---

If you have reached the end of **PART I**, I invite you to continue with the 👉 [Part II – Somatic analysis](README_Part2-3_somatic_analysis.md), where the DNA analysis will be explained step-by-step.

Go back to the top of 👉 [Part I – Preparation & setup](README_Part1-3_setup.md#part-i--preparation--setup)

Go to the main page 👉 [Bash_pipeline_NGS](README.md)
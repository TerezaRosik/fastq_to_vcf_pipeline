{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Genomics Pipeline: FASTQ to VCF\
================================\
\
This repository contains a pipeline for processing genomic data from raw FASTQ files to variant call format (VCF) files. The pipeline supports batch processing, quality control, alignment, duplicate marking, base quality score recalibration (BQSR), and variant calling.\
\
--------------------------------\
Pipeline Overview\
--------------------------------\
\
Scripts:\
1. `fastq_vcf_pipeline.sh`\
   - Main pipeline script for processing multiple samples in batch mode.\
   - Includes optional trimming and BQSR.\
   - Supports SLURM task arrays for parallelisation.\
\
2. `fastq_vcf_pipeline_slurm_job.sh`\
   - SLURM job submission script for running the main pipeline on a cluster.\
\
--------------------------------\
Requirements\
--------------------------------\
\
Software:\
- bwa (for alignment)\
- samtools (for BAM processing)\
- gatk (for variant calling and BQSR)\
- fastqc (for quality control)\
- fastp (for optional trimming)\
- multiqc (for aggregating QC reports)\
- tabix (for indexing VCF files)\
\
Reference Files:\
- Reference genome in FASTA format (e.g., genome.fasta).\
- Known sites VCF files for BQSR (optional).\
\
--------------------------------\
Usage\
--------------------------------\
\
1. Main Pipeline:\
Run the `fastq_vcf_pipeline.sh` script to process multiple samples:\
```\
bash fastq_vcf_pipeline.sh \\\
  -r /path/to/reference.fa \\\
  -s /path/to/samples.tsv \\\
  -o /path/to/output_directory \\\
  -t 4 \\\
  -m 16 \\\
  --trim \\\
  --bqsr \\\
  --known-sites /path/to/known_sites.vcf.gz\
```\
\
Options:\
- `-r`: Path to the reference genome.\
- `-s`: Path to a TSV file with sample information (sample_id, r1, r2).\
- `-o`: Output directory.\
- `-t`: Number of threads.\
- `-m`: Memory in GB.\
- `--trim`: Enable trimming with fastp.\
- `--bqsr`: Enable BQSR with known sites.\
- `--known-sites`: Specify known sites VCF files for BQSR.\
\
2. SLURM Job Submission:\
Submit the pipeline as a SLURM job using the `fastq_vcf_pipeline_slurm_job.sh` script:\
```\
sbatch fastq_vcf_pipeline_slurm_job.sh\
```\
\
SLURM Script Configuration:\
The SLURM script includes the following parameters:\
- `--job-name`: Name of the job.\
- `--output`: Log file for the job.\
- `--cpus-per-task`: Number of CPU cores per task.\
- `--mem`: Memory per node.\
- `--time`: Time limit for the job.\
- `--partition`: Partition name (adjust as needed).\
\
You can modify these parameters in the `fastq_vcf_pipeline_slurm_job.sh` file to suit your cluster environment.\
\
--------------------------------\
Input File Format\
--------------------------------\
\
Samples TSV File:\
The pipeline requires a tab-delimited TSV file with the following columns:\
```\
sample_id    r1    r2\
```\
- `sample_id`: Unique identifier for the sample.\
- `r1`: Path to the forward read FASTQ file.\
- `r2`: Path to the reverse read FASTQ file.\
\
Example:\
```\
sample_id    r1                              r2\
sample1      /path/to/sample1_R1.fastq.gz   /path/to/sample1_R2.fastq.gz\
sample2      /path/to/sample2_R1.fastq.gz   /path/to/sample2_R2.fastq.gz\
```\
\
--------------------------------\
Output Structure\
--------------------------------\
\
The pipeline generates the following directory structure:\
```\
<output_directory>/\
  \uc0\u9500 \u9472 \u9472  logs/               # Log files for each sample\
  \uc0\u9500 \u9472 \u9472  fastqc/             # FastQC reports\
  \uc0\u9500 \u9472 \u9472  trimmed/            # Trimmed FASTQ files (if trimming is enabled)\
  \uc0\u9500 \u9472 \u9472  aligned/            # Sorted BAM files\
  \uc0\u9500 \u9472 \u9472  metrics/            # Alignment and duplication metrics\
  \uc0\u9500 \u9472 \u9472  vcf/                # gVCF files\
  \uc0\u9500 \u9472 \u9472  temp/               # Temporary files\
  \uc0\u9492 \u9472 \u9472  multiqc/            # Aggregated QC reports (if MultiQC is available)\
```\
\
--------------------------------\
Notes\
--------------------------------\
\
- Ensure all required tools are installed and available in your PATH.\
- Use absolute paths for reference files and input directories to avoid errors.\
- Monitor resource usage during the pipeline run to optimise SLURM parameters for larger datasets.\
\
--------------------------------\
Contact\
--------------------------------\
\
For questions or issues, please contact Tereza Rosikova at t.rosikova@uea.ac.uk.}
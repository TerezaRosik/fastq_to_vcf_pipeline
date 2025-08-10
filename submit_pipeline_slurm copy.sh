
#!/bin/bash
#SBATCH --job-name=main_pipeline       # Job name
#SBATCH --output=main_pipeline_%j.log  # Standard output and error log
#SBATCH --ntasks=1                     # Run a single task
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=16G                      # Memory per node
#SBATCH --time=12:00:00                # Time limit (hh:mm:ss)
#SBATCH --partition=standard           # Partition name (adjust as needed)

# Load required modules
module load bwa samtools gatk fastqc fastp multiqc

# Set the path to the main pipeline script
MAIN_PIPELINE_SCRIPT="path/to/pipeline/fastq_vcf_pipeline.sh"

# Run the main pipeline
bash "$MAIN_PIPELINE_SCRIPT" \
  -r /path/to/reference.fa \
  -s /path/to/samples.tsv \
  -o /path/to/output_directory \
  -t 4 \
  -m 16 \
  --trim \
  --bqsr \
  --known-sites /path/to/known_sites.vcf.gz
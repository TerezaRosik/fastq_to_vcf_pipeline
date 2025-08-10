#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# Defaults (override via flags)
THREADS=${THREADS:-4}
MEM_GB=${MEM_GB:-16}
PLATFORM=${PLATFORM:-ILLUMINA}
PLOIDY=${PLOIDY:-2}
DO_TRIM=${DO_TRIM:-0}            # 1 to enable fastp
DO_BQSR=${DO_BQSR:-0}            # 1 to enable BQSR
JAVA_OPTS=${JAVA_OPTS:---java-options "-Xmx${MEM_GB}g"}
TMPDIR_DEFAULT=${TMPDIR_DEFAULT:-"/tmp"}  # Default to system temp directory

# Paths (override via flags)
REFERENCE=""
INPUT_DIR=""
OUTDIR="fastq2gvcf_run_$(date +%Y%m%d_%H%M%S)_$RANDOM"
SAMPLES_TSV=""                   # Tab-delimited: sample_id\tr1\tr2
TASK_ID=""                       # For SLURM array (1-based line index excluding header)
KNOWN_SITES=( )                  # Array of VCF(.gz) for BQSR
ERROR_LOG="$OUTDIR/error.log"    # Centralized error log

usage() {
  cat <<EOF
Usage: $0 \
  -r /path/ref.fa \
  [-i /path/input_fastqs | -s samples.tsv [--task-id N]] \
  [-o /path/outdir] [-t THREADS] [-m MEM_GB] \
  [--trim] [--ploidy P] [--platform ILLUMINA] \
  [--bqsr --known-sites ks1.vcf.gz --known-sites ks2.vcf.gz]

Notes:
  * samples.tsv columns: sample_id\tr1\tr2 (header required)
  * With --task-id (1-based), only that data row is processed (for SLURM arrays).
  * If -i is used without -s, pairs are discovered as *_R1*.f*q.gz <-> *_R2*.f*q.gz.
EOF
}

log(){ echo "[$(date +%F' '%T)] $*" | tee -a "$ERROR_LOG"; }

# Parse arguments
args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference) REFERENCE="$2"; shift 2;;
    -i|--input) INPUT_DIR="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -s|--samples) SAMPLES_TSV="$2"; shift 2;;
    --task-id)    TASK_ID="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -m|--mem)     MEM_GB="$2"; JAVA_OPTS=--java-options\ "-Xmx${MEM_GB}g"; shift 2;;
    --trim)       DO_TRIM=1; shift;;
    --ploidy)     PLOIDY="$2"; shift 2;;
    --platform)   PLATFORM="$2"; shift 2;;
    --bqsr)       DO_BQSR=1; shift;;
    --known-sites) KNOWN_SITES+=("$2"); shift 2;;
    -h|--help) usage; exit 0;;
    *) args+=("$1"); shift;;
  esac
done

[[ -z "$REFERENCE" ]] && { echo "ERROR: --reference is required"; usage; exit 1; }

# Check executables
need() { command -v "$1" >/dev/null 2>&1 || { log "ERROR: $1 not found in PATH"; exit 1; }; }
for exe in bwa samtools gatk fastqc; do need "$exe"; done
if [[ $DO_TRIM -eq 1 ]]; then need fastp; fi
need multiqc || log "WARN: multiqc not in PATH; final aggregation will be skipped"
need tabix || log "WARN: tabix not in PATH; gVCF indexing may fail"

# Prepare output layout
mkdir -p "$OUTDIR"/{logs,fastqc,trimmed,aligned,metrics,vcf,temp,sam}
TMPDIR="${TMPDIR:-$TMPDIR_DEFAULT}"
mkdir -p "$TMPDIR"
touch "$ERROR_LOG"

# Prepare reference indices/dicts
ref_fa="$REFERENCE"
DICT_FILE="${ref_fa%.*}.dict"
if [[ ! -f "${ref_fa}.fai" ]]; then log "Indexing reference (.fai)"; samtools faidx "$ref_fa"; fi
if [[ ! -f "$DICT_FILE" ]]; then log "Creating sequence dictionary"; gatk CreateSequenceDictionary -R "$ref_fa" -O "$DICT_FILE"; fi
if [[ ! -f "${ref_fa}.bwt" ]]; then log "Building BWA index"; bwa index "$ref_fa"; fi

# Helper: discover pairs in a directory
discover_pairs(){
  local dir="$1"
  find "$dir" -type f \( -iname "*R1*.f*q.gz" -o -iname "*_1.f*q.gz" \) | sort |
  while read -r r1; do
    local r2
    r2=$(echo "$r1" | sed -E 's/(_R1_|_1\.)(.*)/\1R2_\2/')
    if [[ -f "$r2" ]]; then
      local base
      base=$(basename "$r1")
      local sample
      sample=${base%%_*}
      echo -e "${sample}\t${r1}\t${r2}"
    else
      log "WARN: Could not find matching R2 for $r1"
    fi
  done
}

# Build an iterable table: sample_id\tr1\tr2
WORK_TSV="$OUTDIR/samples.resolved.tsv"
if [[ -n "$SAMPLES_TSV" ]]; then
  awk 'NR==1{next} {print $1"\t"$2"\t"$3}' "$SAMPLES_TSV" > "$WORK_TSV"
elif [[ -n "$INPUT_DIR" ]]; then
  discover_pairs "$INPUT_DIR" > "$WORK_TSV"
else
  echo "ERROR: Provide either -s samples.tsv or -i input_dir"; exit 1
fi

TOTAL=$(wc -l < "$WORK_TSV" | awk '{print $1}')
if [[ "$TOTAL" -eq 0 ]]; then echo "ERROR: no samples discovered"; exit 1; fi

START=1; END=$TOTAL
if [[ -n "$TASK_ID" ]]; then
  START=$TASK_ID; END=$TASK_ID
  if (( TASK_ID < 1 || TASK_ID > TOTAL )); then echo "ERROR: --task-id $TASK_ID out of range (1..$TOTAL)"; exit 1; fi
fi

# Per-sample processing
process_one(){
  local sample="$1" r1="$2" r2="$3"
  log "=== Processing ${sample} ==="
  local logpref="$OUTDIR/logs/${sample}"

  fastqc -t "$THREADS" -o "$OUTDIR/fastqc" "$r1" "$r2" >"${logpref}.fastqc.log" 2>&1 || { log "FastQC failed"; return 1; }

  local r1_use="$r1" r2_use="$r2"
  if [[ $DO_TRIM -eq 1 ]]; then
    r1_use="$OUTDIR/trimmed/${sample}_R1.trim.fastq.gz"
    r2_use="$OUTDIR/trimmed/${sample}_R2.trim.fastq.gz"
    fastp -i "$r1" -I "$r2" -o "$r1_use" -O "$r2_use" \
          -w "$THREADS" -h "$OUTDIR/trimmed/${sample}.fastp.html" -j "$OUTDIR/trimmed/${sample}.fastp.json" \
          >"${logpref}.fastp.log" 2>&1
  fi

  local sam="$OUTDIR/sam/${sample}.sam"
  local rg="@RG\tID:${sample}\tSM:${sample}\tPL:${PLATFORM}"
  bwa mem -t "$THREADS" -R "$rg" "$ref_fa" "$r1_use" "$r2_use" > "$sam" 2>"${logpref}.bwa.log"

  local bam="$OUTDIR/aligned/${sample}.sorted.bam"
  samtools sort -@ "$THREADS" -o "$bam" "$sam" 2>"${logpref}.sort.log"
  samtools index "$bam"
  samtools flagstat "$bam" > "$OUTDIR/metrics/${sample}.flagstat.txt"
  samtools stats "$bam" > "$OUTDIR/metrics/${sample}.samtools.stats"
  rm -f "$sam"

  local dedup="$OUTDIR/aligned/${sample}.dedup.bam"
  gatk $JAVA_OPTS MarkDuplicatesSpark \
    -I "$bam" -O "$dedup" -M "$OUTDIR/metrics/${sample}.dup_metrics.txt" \
    --tmp-dir "$TMPDIR" --conf 'spark.executor.cores=1' \
    >"${logpref}.markdup.log" 2>&1
  samtools index "$dedup"

  local recal_table="$OUTDIR/metrics/${sample}.bqsr.table"
  local bam_for_call="$dedup"
  if [[ $DO_BQSR -eq 1 && ${#KNOWN_SITES[@]} -gt 0 ]]; then
    gatk $JAVA_OPTS BaseRecalibrator -R "$ref_fa" -I "$dedup" \
      $(printf ' --known-sites %q' "${KNOWN_SITES[@]}") \
      -O "$recal_table" >"${logpref}.bqsr1.log" 2>&1
    local bqsr_bam="$OUTDIR/aligned/${sample}.dedup.bqsr.bam"
    gatk $JAVA_OPTS ApplyBQSR -R "$ref_fa" -I "$dedup" --bqsr-recal-file "$recal_table" -O "$bqsr_bam" \
      >"${logpref}.bqsr2.log" 2>&1
    samtools index "$bqsr_bam"
    bam_for_call="$bqsr_bam"
  fi

  local gvcf="$OUTDIR/vcf/${sample}.g.vcf.gz"
  gatk $JAVA_OPTS HaplotypeCaller -R "$ref_fa" -I "$bam_for_call" -O "$gvcf" -ERC GVCF \
    --native-pair-hmm-threads "$THREADS" --sample-ploidy "$PLOIDY" \
    --tmp-dir "$TMPDIR" >"${logpref}.hc.log" 2>&1

  tabix -f "$gvcf" 2>/dev/null || true
  log "=== Completed ${sample} ==="
}

ROW=0
while IFS=$'\t' read -r sample r1 r2; do
  ROW=$((ROW+1))
  if (( ROW < START || ROW > END )); then continue; fi
  process_one "$sample" "$r1" "$r2"
  sync || true
done < "$WORK_TSV"

if command -v multiqc >/dev/null 2>&1; then
  (cd "$OUTDIR" && multiqc . -o "$OUTDIR/multiqc") || log "WARN: MultiQC failed"
fi

log "All done. Outputs in $OUTDIR"
#!/usr/bin/env bash
# ==============================================================================
# 03_primer_trimming.sh
# Removes 16S V4 primers (515F / 806R, Earth Microbiome Project primers)
# from paired-end reads with cutadapt. Read pairs without a detectable primer
# (i.e. off-target amplification) are discarded.
# ==============================================================================
set -euo pipefail

RAW_DIR="$(dirname "$0")/../data/raw"
TRIM_DIR="$(dirname "$0")/../results/otu/trimmed"
LOG_DIR="$(dirname "$0")/../results/qc/cutadapt_logs"
mkdir -p "${TRIM_DIR}" "${LOG_DIR}"

FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"   # 515F
REV_PRIMER="GGACTACNVGGGTWTCTAAT"  # 806R

for r1 in "${RAW_DIR}"/*_R1.fastq.gz; do
  sample=$(basename "${r1}" _R1.fastq.gz)
  r2="${RAW_DIR}/${sample}_R2.fastq.gz"

  echo ">> Trimming primers for ${sample}"
  cutadapt \
    -g "${FWD_PRIMER}" \
    -G "${REV_PRIMER}" \
    --discard-untrimmed \
    --minimum-length 50 \
    -o "${TRIM_DIR}/${sample}_R1.trimmed.fastq.gz" \
    -p "${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz" \
    "${r1}" "${r2}" \
    > "${LOG_DIR}/${sample}.cutadapt.log" 2>&1
done

echo ">> Primer trimming complete. Summary (reads retained):"
for log in "${LOG_DIR}"/*.log; do
  sample=$(basename "${log}" .cutadapt.log)
  pairs=$(grep "Pairs written" "${log}" | head -1)
  echo "   ${sample}: ${pairs}"
done

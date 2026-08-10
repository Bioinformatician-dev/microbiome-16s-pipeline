#!/usr/bin/env bash
# ==============================================================================
# 02_quality_control.sh
# Runs FastQC on all raw FASTQ files and summarizes with MultiQC.
# ==============================================================================
set -euo pipefail

RAW_DIR="$(dirname "$0")/../data/raw"
QC_DIR="$(dirname "$0")/../results/qc"
mkdir -p "${QC_DIR}/fastqc_raw"

echo ">> Running FastQC on raw reads..."
fastqc -q -o "${QC_DIR}/fastqc_raw" "${RAW_DIR}"/*.fastq.gz

echo ">> Summarizing with MultiQC..."
multiqc -q -f -o "${QC_DIR}" "${QC_DIR}/fastqc_raw" >/dev/null 2>&1 || \
  python3 -m multiqc -q -f -o "${QC_DIR}" "${QC_DIR}/fastqc_raw"

echo ">> QC complete. Reports in ${QC_DIR}"

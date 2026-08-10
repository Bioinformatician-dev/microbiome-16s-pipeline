#!/usr/bin/env bash
# ==============================================================================
# 04_merge_filter.sh
# Merges primer-trimmed paired-end reads into full-length V4 amplicons and
# applies quality filtering (max expected error, no ambiguous bases),
# using VSEARCH (Rognes et al. 2016, PeerJ, doi:10.7717/peerj.2584).
# ==============================================================================
set -euo pipefail

TRIM_DIR="$(dirname "$0")/../results/otu/trimmed"
MERGE_DIR="$(dirname "$0")/../results/otu/merged"
LOG_DIR="$(dirname "$0")/../results/qc/vsearch_logs"
mkdir -p "${MERGE_DIR}" "${LOG_DIR}"

for r1 in "${TRIM_DIR}"/*_R1.trimmed.fastq.gz; do
  sample=$(basename "${r1}" _R1.trimmed.fastq.gz)
  r2="${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz"

  echo ">> Merging pairs for ${sample}"
  vsearch --fastq_mergepairs "${r1}" --reverse "${r2}" \
    --fastqout "${MERGE_DIR}/${sample}.merged.fastq" \
    --fastq_minovlen 20 \
    --fastq_maxdiffpct 15 \
    --relabel "${sample}." \
    --threads 2 \
    2> "${LOG_DIR}/${sample}.merge.log"

  echo ">> Quality-filtering merged reads for ${sample}"
  vsearch --fastq_filter "${MERGE_DIR}/${sample}.merged.fastq" \
    --fastq_maxee 1.0 \
    --fastq_minlen 200 \
    --fastq_maxlen 300 \
    --fastq_maxns 0 \
    --fastaout "${MERGE_DIR}/${sample}.filtered.fasta" \
    --threads 2 \
    2>> "${LOG_DIR}/${sample}.merge.log"
done

echo ">> Merge/filter complete. Read survival per sample:"
for log in "${LOG_DIR}"/*.merge.log; do
  sample=$(basename "${log}" .merge.log)
  merged=$(grep -m1 "Merged" "${log}" || true)
  kept=$(grep -m1 "sequences kept" "${log}" || true)
  echo "   ${sample}: ${merged} | ${kept}"
done

echo ">> Pooling all samples into one FASTA for dereplication..."
cat "${MERGE_DIR}"/*.filtered.fasta > "${MERGE_DIR}/all_samples.filtered.fasta"
grep -c ">" "${MERGE_DIR}/all_samples.filtered.fasta" | xargs echo "   total filtered reads pooled:"

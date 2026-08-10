#!/usr/bin/env bash
# ==============================================================================
# 05_denoise_cluster.sh
# Dereplicates pooled reads, denoises them into amplicon sequence variants
# (ASVs / zOTUs) with the UNOISE3 algorithm, removes chimeras de novo, and
# builds a sample-by-ASV abundance (OTU) table.
# ==============================================================================
set -euo pipefail

MERGE_DIR="$(dirname "$0")/../results/otu/merged"
OTU_DIR="$(dirname "$0")/../results/otu"
mkdir -p "${OTU_DIR}"

ALL_FASTA="${MERGE_DIR}/all_samples.filtered.fasta"

echo ">> Dereplicating pooled sequences..."
vsearch --derep_fulllength "${ALL_FASTA}" \
  --sizeout \
  --relabel Uniq. \
  --output "${OTU_DIR}/uniques.fasta" \
  --threads 2 2> "${OTU_DIR}/derep.log"

echo ">> Denoising with UNOISE3 (ASV inference, min size 4)..."
vsearch --cluster_unoise "${OTU_DIR}/uniques.fasta" \
  --minsize 4 \
  --unoise_alph 2 \
  --centroids "${OTU_DIR}/asvs_raw.fasta" \
  --threads 2 2> "${OTU_DIR}/denoise.log"

echo ">> Removing chimeras de novo (UCHIME3)..."
vsearch --uchime3_denovo "${OTU_DIR}/asvs_raw.fasta" \
  --nonchimeras "${OTU_DIR}/asvs.fasta" \
  --relabel ASV_ \
  2> "${OTU_DIR}/chimera.log"

n_asv=$(grep -c ">" "${OTU_DIR}/asvs.fasta")
echo ">> ${n_asv} chimera-free ASVs retained."

echo ">> Mapping reads from each sample back to ASVs to build abundance table..."
vsearch --usearch_global "${ALL_FASTA}" \
  --db "${OTU_DIR}/asvs.fasta" \
  --id 0.97 \
  --otutabout "${OTU_DIR}/otu_table.tsv" \
  --threads 2 2> "${OTU_DIR}/otutab.log"

echo ">> Final ASV/OTU table: ${OTU_DIR}/otu_table.tsv"
head -1 "${OTU_DIR}/otu_table.tsv"
wc -l "${OTU_DIR}/otu_table.tsv"

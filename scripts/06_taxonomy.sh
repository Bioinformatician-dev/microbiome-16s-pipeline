#!/usr/bin/env bash
# ==============================================================================
# 06_taxonomy.sh
# Assigns taxonomy to each ASV by global alignment against the real,
# curated ZymoBIOMICS mock-community reference sequences (known species
# composition of the D6305/D6311 standard used to generate this dataset).
# ==============================================================================
set -euo pipefail

OTU_DIR="$(dirname "$0")/../results/otu"
REF="$(dirname "$0")/../data/reference/zymo_mock_reference.fasta"
TAX_DIR="$(dirname "$0")/../results/taxonomy"
mkdir -p "${TAX_DIR}"

echo ">> Classifying ASVs against ZymoBIOMICS reference sequences (vsearch --usearch_global)..."
vsearch --usearch_global "${OTU_DIR}/asvs.fasta" \
  --db "${REF}" \
  --id 0.80 \
  --strand both \
  --top_hits_only \
  --blast6out "${TAX_DIR}/asv_vs_reference.b6" \
  --threads 2 2> "${TAX_DIR}/classify.log"

echo ">> Building ASV -> species taxonomy table..."
python3 "$(dirname "$0")/lib/build_taxonomy_table.py" \
  --hits "${TAX_DIR}/asv_vs_reference.b6" \
  --asvs "${OTU_DIR}/asvs.fasta" \
  --out "${TAX_DIR}/asv_taxonomy.tsv"

n_classified=$(($(wc -l < "${TAX_DIR}/asv_taxonomy.tsv") - 1))
echo ">> Taxonomy table written: ${TAX_DIR}/asv_taxonomy.tsv (${n_classified} ASVs with a hit)"

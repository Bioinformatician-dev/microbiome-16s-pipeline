#!/usr/bin/env bash
# ==============================================================================
# run_pipeline.sh — runs the full 16S rRNA microbial-community pipeline
# end to end: download -> QC -> primer trim -> merge/filter -> denoise/OTU ->
# taxonomy -> diversity analysis.
# ==============================================================================
set -euo pipefail
cd "$(dirname "$0")"

echo "############################################"
echo "# 1/7 Downloading real 16S amplicon data"
echo "############################################"
bash scripts/01_download_data.sh

echo "############################################"
echo "# 2/7 Raw read quality control"
echo "############################################"
bash scripts/02_quality_control.sh

echo "############################################"
echo "# 3/7 Primer trimming"
echo "############################################"
bash scripts/03_primer_trimming.sh

echo "############################################"
echo "# 4/7 Paired-end merging & quality filtering"
echo "############################################"
bash scripts/04_merge_filter.sh

echo "############################################"
echo "# 5/7 Denoising (ASVs) & chimera removal"
echo "############################################"
bash scripts/05_denoise_cluster.sh

echo "############################################"
echo "# 6/7 Taxonomic classification"
echo "############################################"
bash scripts/06_taxonomy.sh

echo "############################################"
echo "# 7/7 Diversity analysis & figures"
echo "############################################"
python3 scripts/07_diversity_analysis.py

echo
echo "Pipeline complete. See results/ for tables and results/figures/ for plots."

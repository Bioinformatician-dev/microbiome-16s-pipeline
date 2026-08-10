#!/usr/bin/env bash
# ==============================================================================
# 01_download_data.sh
#
# Downloads real paired-end Illumina MiSeq 16S rRNA (V4 region, 515F/806R)
# amplicon sequencing data of the ZymoBIOMICS Microbial Community Standard
# (a defined mock community of 8 bacteria + 2 fungi with known composition).
#
# Source: nf-core/test-datasets (ampliseq branch) - real sequencing reads
# used as the reference test data for the nf-core/ampliseq pipeline
# (Straub et al. 2020, Front. Microbiol., doi:10.3389/fmicb.2020.550420).
# The reads derive from real ENA-deposited biological samples
# (materialSampleID accessions ERS3508481-ERS3508484).
#
# Two mock-community abundance layouts are included, each with 2 replicates:
#   - "b" (Std, even distribution)   -> samples 1, 2
#   - "a" (Log, staggered distribution) -> samples 1a, 2a
# ==============================================================================
set -euo pipefail

RAW_DIR="$(dirname "$0")/../data/raw"
REF_DIR="$(dirname "$0")/../data/reference"
META_DIR="$(dirname "$0")/../data/metadata"
mkdir -p "$RAW_DIR" "$REF_DIR" "$META_DIR"

BASE_URL="https://raw.githubusercontent.com/nf-core/test-datasets/ampliseq"

declare -A SAMPLES=(
  [sampleID_1]="testdata/1_S103_L001"
  [sampleID_1a]="testdata/1a_S103_L001"
  [sampleID_2]="testdata/2_S115_L001"
  [sampleID_2a]="testdata/2a_S115_L001"
)

echo ">> Downloading raw paired-end FASTQ files..."
for sample in "${!SAMPLES[@]}"; do
  prefix="${SAMPLES[$sample]}"
  for r in R1 R2; do
    out="${RAW_DIR}/${sample}_${r}.fastq.gz"
    url="${BASE_URL}/${prefix}_${r}_001.fastq.gz"
    echo "   ${sample} ${r}: ${url}"
    curl -sSL --fail "${url}" -o "${out}"
  done
done

echo ">> Downloading ZymoBIOMICS mock-community reference 16S/18S sequences (ground truth)..."
curl -sSL --fail \
  "${BASE_URL}/expected/D6305-D6311_ZymoBIOMICS.STD.refseq.v2.ssrRNAs.unique.fasta" \
  -o "${REF_DIR}/zymo_mock_reference.fasta"

echo ">> Downloading expected (theoretical) mock-community abundances..."
curl -sSL --fail \
  "${BASE_URL}/expected/D6305-D6311_ZymoBIOMICS.STD.refseq.v2.ssrRNAs.unique.abundance.tsv" \
  -o "${REF_DIR}/zymo_expected_abundance.tsv"

echo ">> Downloading sample metadata..."
curl -sSL --fail "${BASE_URL}/samplesheets/Metadata.tsv" -o "${META_DIR}/metadata.tsv"

echo ">> Done. Raw data:"
ls -lh "${RAW_DIR}"
echo ">> Verifying read counts per file..."
for f in "${RAW_DIR}"/*.fastq.gz; do
  n=$(( $(zcat "$f" | wc -l) / 4 ))
  printf "   %-35s %6d reads\n" "$(basename "$f")" "$n"
done

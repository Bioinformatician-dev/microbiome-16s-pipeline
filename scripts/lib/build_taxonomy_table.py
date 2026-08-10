#!/usr/bin/env python3
"""Parse vsearch --blast6out hits against the ZymoBIOMICS reference into a
simple ASV -> genus/species/marker taxonomy table."""
import argparse
import csv


def parse_ref_header(ref_id: str):
    """'Bacillus_subtilis_16S_4x_1-8-9-10' -> genus, species, marker_gene."""
    parts = ref_id.split("_")
    genus, species = parts[0], parts[1]
    marker = "16S" if "16S" in ref_id else ("18S" if "18S" in ref_id else "rRNA")
    return genus, species, marker


def load_asv_ids(fasta_path):
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="vsearch blast6out file")
    ap.add_argument("--asvs", required=True, help="ASV fasta (for full ASV list)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    all_asvs = load_asv_ids(args.asvs)

    best_hit = {}
    with open(args.hits) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            query, target, pident = row[0], row[1], float(row[2])
            # blast6out is already top-hit-only per query (--top_hits_only),
            # keep the first (best) hit seen for each query.
            if query not in best_hit:
                best_hit[query] = (target, pident)

    with open(args.out, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["ASV_ID", "genus", "species", "marker_gene",
                          "reference_hit", "percent_identity"])
        for asv in all_asvs:
            if asv in best_hit:
                target, pident = best_hit[asv]
                genus, species, marker = parse_ref_header(target)
                writer.writerow([asv, genus, species, marker, target, f"{pident:.2f}"])
            else:
                writer.writerow([asv, "Unclassified", "Unclassified", "NA", "NA", "NA"])


if __name__ == "__main__":
    main()

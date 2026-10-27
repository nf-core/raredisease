#!/usr/bin/env python3
"""Parse GATK CalculateContamination output for MultiQC."""

import argparse


def parse_contamination(contamination_table, sample_id, prefix):
    """Parse contamination table and write MultiQC output."""
    with open(contamination_table, 'r') as f:
        lines = f.readlines()
        data_line = lines[1].strip().split('\t')
        contamination = float(data_line[1])
        contamination_pct = contamination * 100

    with open(f"{prefix}_contamination_mqc.tsv", 'w') as out:
        out.write("# id: 'gatk_contamination'\n")
        out.write("# section_name: 'GATK Contamination'\n")
        out.write("# description: 'Sample contamination estimates from GATK CalculateContamination'\n")
        out.write("# plot_type: 'generalstats'\n")
        out.write("# pconfig:\n")
        out.write("#     contamination_pct:\n")
        out.write("#         title: 'Contamination'\n")
        out.write("#         description: 'Estimated sample contamination percentage'\n")
        out.write("#         max: 10\n")
        out.write("#         min: 0\n")
        out.write("#         scale: 'RdYlGn-rev'\n")
        out.write("#         suffix: '%'\n")
        out.write("#         format: '{:,.2f}'\n")
        out.write("Sample\tcontamination_pct\n")
        out.write(f"{sample_id}\t{contamination_pct:.4f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse GATK contamination output for MultiQC")
    parser.add_argument("--input", required=True, help="GATK contamination table")
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    args = parser.parse_args()

    parse_contamination(args.input, args.sample_id, args.prefix)

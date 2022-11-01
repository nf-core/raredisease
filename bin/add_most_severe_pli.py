#!/usr/bin/env python
import argparse
import sys
from pathlib import Path
from typing import TextIO


def parse_vep_transcripts(transcripts: list, symbol_ind: int) -> list:
    """
    Parse each transcript and return a list of gene symbols.

    Args:
        transcripts (list): A list of vep transcript annotation
        symbol_ind  (int) : Index of the "allele" in the vep annotation record

    Returns:
        gene_ids     (list): list of gene ids in the record
    """

    gene_ids = []
    for transcript in transcripts:
        vep_fields = transcript.strip().split("|")
        gene_id = vep_fields[symbol_ind]
        gene_ids.append(gene_id)
    return gene_ids


def construct_most_severe_pli_info(line: str, symbol_ind: int, pli_gene: dict) -> list:
    """
    Parse gene symbols, find the highest pli value of all gene symbols, add most_severe_pli tag to the info
    field and return a list of modified columns

    Args:
        line        (str) : Vcf record
        symbol_ind  (int) : Index of the "SYMBOL" in the vep annotation record
        pli_gene    (dict): A dict of pli values, where gene symbols are the keys

    Returns:
        columns     (list): A list of fields in the vcf record with most severe pli added
                            to the INFO column
    """

    columns = line.strip().split()
    info_fields = columns[7].split(";")
    for field in info_fields:
        if field.startswith("CSQ="):
            transcripts = field.split("CSQ=")[1].split(",")
    gene_ids = parse_vep_transcripts(transcripts, symbol_ind)
    unique_ids = list(set(gene_ids))
    pli_values = []
    for gene_id in unique_ids:
        if gene_id != "" and pli_gene.get(gene_id) is not None:
            pli_values.append(pli_gene.get(gene_id))
    if pli_values:
        columns[7] += ";most_severe_pli={:.2f}".format(max(pli_values))
    return columns


def parse_vep_csq_schema(line: str) -> int:
    """
    Get indices of gene symbol in the annotation

    Args:
        line: INFO line in the vcf header with CSQ information

    Returns:
        symbol_ind  (int) : Index of the "SYMBOL" in the vep annotation record
    """
    fields = line.strip().split("Format: ")[1].replace('">', "").split("|")
    symbol_ind = fields.index("SYMBOL")

    return symbol_ind


def write_pli_annotated_vcf(file_in: TextIO, file_out: TextIO, var_csq: list):
    """Add most severe pli field to record, and write the record to a vcf file"""
    for line in file_in:
        if line.startswith("#"):
            file_out.write(line)
            if line.startswith("##INFO=<ID=CSQ"):
                symbol_ind = parse_vep_csq_schema(line)
                file_out.write(
                    '##INFO=<ID=most_severe_pli,Number=1,Type=Float,Description="Probabililty of a gene being loss-of-function intolerant score.">\n'
                )
        else:
            vcf_record = construct_most_severe_pli_info(line, symbol_ind, var_csq)
            file_out.write("\t".join(vcf_record) + "\n")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate vcf with the most severe pli field.",
        epilog="Example: python vcfparser.py --file_in vep.vcf --file_out vep.most_severe_pli.vcf --pli pli_per_gene.txt",
    )
    parser.add_argument(
        "--file_in",
        metavar="FILE_IN",
        type=Path,
        help="Vcf file annotated with vep.",
    )
    parser.add_argument(
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Vcf with most_severe_pli annotations added to it.",
    )
    parser.add_argument(
        "--pli",
        metavar="PLI",
        type=Path,
        help="Pli",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if not args.file_in.is_file():
        print(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    if not args.pli.is_file():
        print(f"The given variant consequence file {args.pli} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    pli_gene = {}
    with open(args.pli) as f:
        for line in f:
            cols = line.strip().split()
            if cols[0] != "gene":
                pli_gene[cols[0]] = float(cols[1])
    with open(args.file_out, "w") as out_vcf:
        with open(args.file_in) as in_vcf:
            write_pli_annotated_vcf(in_vcf, out_vcf, pli_gene)


if __name__ == "__main__":
    sys.exit(main())

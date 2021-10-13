#!/usr/bin/env python

import os
import sys
import gzip
import argparse

def parse_args(args=None):
    Description = "Check that input vcf files have been normalized."
    Epilog = "Example usage: python check_input_vcf.py --INPUT_VCFS <VCF_FILE_1> <VCF_FILE_2>...<VCF_FILE_N> --OUTPUT <output.txt>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--INPUT_VCFS", help="Input vcf files.", nargs= '*')
    parser.add_argument("--OUTPUT", help="Output file containing a list of files that need to be normalized.")
    return parser.parse_args(args)


def check_vcf(files_in, file_out):
    """
    This function checks that input vcf files have been normalized with bcftools, and writes the filenames of
    those that need to be normalized to a text file

    """

    with open(file_out,'w') as out:
        for file in files_in:
            if file.endswith(".gz"):
                with gzip.open(file,'rt') as vcf:
                    for line in vcf:
                        if line.startswith("##bcftools_norm"):
                            break
                        elif not line.startswith("#"):
                            base = os.path.basename(file).split(".")[0]
                            out.write(base + "," + os.path.abspath(file) + "\n")
                            break
            else:
                print("Please compress %s using bgzip" %file )




def main(args=None):
    args = parse_args(args)
    check_vcf(args.INPUT_VCFS, args.OUTPUT)


if __name__ == "__main__":
    sys.exit(main())

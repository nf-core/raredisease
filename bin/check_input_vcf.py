#!/usr/bin/env python

import os
import sys
import gzip
import argparse

def parse_args(args=None):
    Description = "Check that input vcf file has been normalized."
    Epilog = "Example usage: python check_input_vcf.py --INPUT_VCF <VCF_FILE_1> <VCF_FILE_2>...<VCF_FILE_N> --OUTPUT <output.txt>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--INPUT_VCF", help="Input vcf file.")
    parser.add_argument("--OUTPUT", help="Output file containing a list of files that need to be normalized.")
    return parser.parse_args(args)


def check_vcf(file_in, file_out):
    """
    This function checks that an input vcf file has been normalized with bcftools, and if not, it writes its basename and the filename
    to a text file

    """

    with open(file_out,'w') as out:
        if file_in.endswith(".gz"):
            with gzip.open(file_in,'rt') as vcf:
                for line in vcf:
                    if line.startswith("##bcftools_norm"):
                        base = os.path.basename(file_in).rsplit(".",2)[0]
                        out.write("id,filepath,processed\n")
                        out.write(base + "," + os.path.abspath(file_in) + ",yes\n")
                        break
                    elif not line.startswith("#"):
                        base = os.path.basename(file_in).rsplit(".",2)[0]
                        out.write("id,filepath,processed\n")
                        out.write(base + "," + os.path.abspath(file_in) + ",no\n")
                        break
        else:
            print("Please compress %s using bgzip" %file_in )


def main(args=None):
    args = parse_args(args)
    check_vcf(args.INPUT_VCF, args.OUTPUT)


if __name__ == "__main__":
    sys.exit(main())

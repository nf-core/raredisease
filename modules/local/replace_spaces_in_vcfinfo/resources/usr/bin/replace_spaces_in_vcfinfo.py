#!/usr/bin/env python3


# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Reformat VCF by replacing spaces with underscores in INFO field"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input VCF file"
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_fn = input_path.stem + "_reformatted.vcf"

    with input_path.open() as input_file, open(output_fn, "w") as output:
        for line in input_file:
            if line.startswith("#"):
                output.write(line)
            else:
                spl = line.rstrip("\n").split("\t")
                if len(spl) > 7:
                    spl[7] = spl[7].replace(" ", "_")
                output.write("\t".join(spl) + "\n")


if __name__ == "__main__":
    main()

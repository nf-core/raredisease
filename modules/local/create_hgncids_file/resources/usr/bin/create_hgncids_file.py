#!/usr/bin/env python3

# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Reformat input file based on meta.id"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input file"
    )
    parser.add_argument(
        "--meta-id",
        required=True,
        help="Value of meta.id (e.g. scout)"
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_fn = input_path.stem + "_reformatted.txt"

    with input_path.open() as input_file, open(output_fn, "w") as output:
        if args.meta_id == "scout":
            for line in input_file:
                if not line.startswith("#") and line.strip():
                    spl = line.rstrip("\n").split("\t")
                    if len(spl) > 3:
                        output.write(spl[3] + "\n")
        else:
            for line in input_file:
                if not line.startswith("#"):
                    output.write(line)

if __name__ == "__main__":
    main()

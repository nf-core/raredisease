#!/usr/bin/env python3

# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

import argparse
import gzip
import re
import sys
from pathlib import Path


INDEX_SUFFIXES = {".tbi", ".csi", ".bai"}


def has_records(path: Path) -> bool:
    if path.suffix in INDEX_SUFFIXES:
        return True
    open_fn = gzip.open if path.suffix in {".gz", ".bgz"} else open
    with open_fn(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                return True
    return False


def split_toml(content: str):
    """Return (header, blocks) where blocks is a list of [[annotation]] block strings."""
    parts = re.split(r"(?=\[\[annotation\]\])", content)
    return parts[0], parts[1:]


def block_filename(block: str) -> str | None:
    match = re.search(r'^file\s*=\s*["\']([^"\']+)["\']', block, re.MULTILINE)
    if match:
        return Path(match.group(1)).name
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Check vcfanno database files and remove zero-record entries from the TOML"
    )
    parser.add_argument("--toml", required=True, help="Path to vcfanno TOML config file")
    parser.add_argument(
        "--databases",
        nargs="+",
        required=True,
        help="Paths to vcfanno database files",
    )
    args = parser.parse_args()

    empty_basenames = set()
    for db in args.databases:
        path = Path(db)
        if not has_records(path):
            print(f"WARNING: {path.name} has zero records — removing from TOML", file=sys.stderr)
            empty_basenames.add(path.name)
        else:
            print(f"{path.name}: has records")

    toml_path = Path(args.toml)
    content = toml_path.read_text()
    header, blocks = split_toml(content)

    filtered_blocks = [b for b in blocks if block_filename(b) not in empty_basenames]
    removed = len(blocks) - len(filtered_blocks)
    if removed:
        print(f"Removed {removed} annotation block(s) from TOML.", file=sys.stderr)

    output_path = Path(toml_path.stem + "_filtered.toml")
    output_path.write_text(header + "".join(filtered_blocks))
    print(f"Filtered TOML written to {output_path}")


if __name__ == "__main__":
    main()

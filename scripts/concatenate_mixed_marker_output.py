#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2026 nylander <johan.nylander@nrm.se>
#
# Distributed under terms of the MIT license.

"""
Concatenate marker-specific fasta files into one fasta file per marker. Note:
fasta headers need to have labels following a specific format. For example:
[ref:PV424158.1|ITS|1] "
"""

import argparse
import textwrap
import os
import re
from collections import defaultdict


DEFAULT_EXCLUDED_SUFFIXES = (
    ".bcftools.fasta",
    ".samtools-a.fasta",
    ".samtools.fasta",
    ".samtools-iupac.fasta",
    ".samtools-a-iupac.fasta",
)


def parse_fasta_files(input_folder, suffix=None):
    records = defaultdict(list)
    ref_pattern = re.compile(r"\[ref:(.+?)\]")
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta" if suffix is None else suffix) and not (
            suffix is None and filename.endswith(DEFAULT_EXCLUDED_SUFFIXES)
        ):
            with open(os.path.join(input_folder, filename), "r") as file:
                current_key = None
                current_header = None
                current_sequence = []
                for line in file:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if current_key and current_header and current_sequence:
                            records[current_key].append(
                                (current_header, "".join(current_sequence))
                            )
                        current_sequence = []
                        current_header = line[1:]
                        match = ref_pattern.search(line)
                        if match:
                            current_key = match.group(1)
                            if not re.match(r"^[\w\-.]+$", current_key):
                                current_key = re.sub(r"[^\w\-.]", "_", current_key)
                        else:
                            current_key = None
                            current_header = None
                    else:
                        if current_key and current_header:
                            current_sequence.append(line)
                if current_key and current_header and current_sequence:
                    records[current_key].append((current_header, "".join(current_sequence)))
    return records


def write_output_files(records, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for key, rec_list in records.items():
        output_file = os.path.join(output_folder, f"{key}.fasta")
        with open(output_file, "w") as file:
            for header, seq in rec_list:
                file.write(f">{header}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''
            Concatenate marker-specific fasta files into one fasta file per marker.
            Fasta headers need to have labels following a specific format.
            For example: ">... [ref:PV424158.1|ITS|1]".
            By default, reads .fasta files except those ending in .bcftools.fasta,
            .samtools-a.fasta, .samtools.fasta, .samtools-iupac.fasta, or
            .samtools-a-iupac.fasta. An explicit --suffix disables these exclusions.
            '''
        ),formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_folder",
                        help="Path to the input folder containing fasta files")
    parser.add_argument("output_folder",
                        help="Path to the output folder where fasta files will be saved")
    parser.add_argument("--suffix",
                        help="Read all files ending in this suffix (default: .fasta with exclusions)")
    args = parser.parse_args()
    records = parse_fasta_files(args.input_folder, args.suffix)
    write_output_files(records, args.output_folder)

if __name__ == "__main__":
    main()

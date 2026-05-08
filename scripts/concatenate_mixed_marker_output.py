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
import os
import re
from collections import defaultdict

def parse_fasta_files(input_folder):
    records = defaultdict(list)
    ref_pattern = re.compile(r"\[ref:(.+?)\]")
    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta") and not any(
            filename.endswith(suffix)
            for suffix in [
                ".bcftools.fasta",
                ".samtools-a.fasta",
                ".samtools.fasta",
                ".samtools-iupac.fasta",
            ]
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
        description=(
            "Concatenate marker-specific fasta files into one fasta file per marker.\n"
            "Note: fasta headers need to have labels following a specific format.\n"
            "For example: >... [ref:PV424158.1|ITS|1]"
        )
    )
    parser.add_argument("input_folder",
                        help="Path to the input folder containing fasta files")
    parser.add_argument("output_folder",
                        help="Path to the output folder where fasta files will be saved")
    args = parser.parse_args()
    records = parse_fasta_files(args.input_folder)
    write_output_files(records, args.output_folder)

if __name__ == "__main__":
    main()

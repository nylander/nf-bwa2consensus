#! /usr/bin/env python3
# vim:fenc=utf-8

# Copyright © 2023 nylander <johan.nylander@nrm.se>
# Distributed under terms of the MIT license.
# Last modified: 2026-04-10 19:00:42
# Sign: JN

"""
Map fastq sequences to reference using bwa mem, calculate consensus sequence
using samtools consensus and/or bcftools.

Input: file_R1.fq.gz, file_R2.fq.gz, ref.fasta
Output: consensus.fasta

Dependencies: bwa, samtools, bcftools

Samtools consensus documentation: <https://www.htslib.org/doc/1.15/samtools-consensus.html>
Bcftools documentation: <https://samtools.github.io/bcftools/bcftools.html>

Note: The use of bcftools mpileup assumes diploid individuals.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
import logging
from typing import List
#from typing import List, Tuple # use Tuple if several return values: -> Tuple[str, str]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def parse_args(args: List[str]) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Map fastq sequences to reference using bwa mem and calculate "
            "consensus sequence using samtools consensus or bcftools"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-1", "--fastq1",
        type=str,
        required=True,
        help="Fastq file with forward reads (can be compressed)",
    )
    parser.add_argument(
        "-2", "--fastq2",
        type=str,
        required=True,
        help="Fastq file with reverse reads (can be compressed)",
    )
    parser.add_argument(
        "-r", "--ref",
        type=str,
        required=True,
        help="Reference fasta file",
    )
    parser.add_argument(
        "-o", "--outdir",
        type=str,
        default=".",
        help="Output directory",
    )
    parser.add_argument(
        "-p", "--prefix",
        type=str,
        default="consensus",
        help="Prefix for output files",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of threads to use",
    )
    parser.add_argument(
        "--maxdepth",
        type=int,
        default=50,
        help="Maxdepth ...",
    )
    parser.add_argument(
        "--mindepth",
        type=int,
        default=1,
        help="Mindepth ...",
    )
    #parser.add_argument(
    #    "--ploidy",
    #    type=int,
    #    default=2,
    #    help="Ploidy level",
    #)
    return parser.parse_args(args)

def run_bwa_index_of_ref(ref:str) -> None:
    """Run bwa index."""
    cmd = f"bwa index {ref} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)

def run_samtools_index(bam:str, threads:int) -> None:
    """Run samtools index."""
    cmd = f"samtools index -@ {threads} {bam}"
    subprocess.run(cmd, shell=True, check=True)

def run_bwa(fastq1:str, fastq2:str, ref:str, outdir:str, prefix:str, threads:int) -> None:
    """Run bwa mem."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outbam = os.path.join(outdir, f"{prefix}.bam")
    cmd = (
        f"bwa mem -t {threads} {ref} {fastq1} {fastq2} 2>/dev/null"
        f"| samtools sort --threads {threads} -o {outbam} 2>/dev/null"
    )
    subprocess.run(cmd, shell=True, check=True)

def run_samtools_consensus(bam:str, outdir:str, prefix:str, threads:str) -> str:
    """Run samtools consensus."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfasta = os.path.join(outdir, f"{prefix}.samtools.fasta")
    cmd = f"samtools consensus --threads {threads} {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_consensus_a(bam:str, outdir:str, prefix:str, threads:str) -> str:
    """
    Run samtools consensus with -a option (Outputs all bases, from start to end of reference,
    even when the aligned data does not extend to the ends).
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfasta = os.path.join(outdir, f"{prefix}.samtools_a.fasta")
    cmd = f"samtools consensus --threads {threads} -a {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_consensus_iupac(bam:str, outdir:str, prefix:str, threads:str) -> str:
    """Run samtools consensus with --ambig option."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfasta = os.path.join(outdir, f"{prefix}.samtools_iupac.fasta")
    cmd = f"samtools consensus --threads {threads} --ambig {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_stats(bam:str, outdir:str, prefix:str, threads:str) -> None:
    """Run samtools stats."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outstats = os.path.join(outdir, f"{prefix}.stats")
    cmd = f"samtools stats --threads {threads} {bam} > {outstats}"
    subprocess.run(cmd, shell=True, check=True)

def run_samtools_depth_a(bam:str, outdir:str, prefix:str, threads:str) -> None:
    """Run samtools depth with -a option"""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdepth = os.path.join(outdir, f"{prefix}.depth.tsv")
    cmd = f"samtools depth -a --threads {threads} {bam} > {outdepth}"
    subprocess.run(cmd, shell=True, check=True)

def run_bcftools_mpileup(ref:str, maxdepth:str, bam:str,  mindepth:str, tmp_vcf_gz:str) -> None:
    """Run bcftools mpileup"""
    cmd = (
        f"bcftools mpileup --fasta-ref {ref} --annotate DP "
        f"--max-depth {maxdepth} {bam} "
        "| bcftools call --multiallelic-caller --output-type z "
        "| bcftools view --max-alleles 2 "
        f'--include " INFO/INDEL=0 && FORMAT/DP>={mindepth} " '
        f"--output-type z --output {tmp_vcf_gz}"
    )
    subprocess.run(cmd, shell=True, check=True)

def run_bcftools_index(tmp_vcf_gz:str, threads:int) -> None:
    """Run bcftools index"""
    cmd = f"bcftools index --threads {threads} {tmp_vcf_gz}"
    subprocess.run(cmd, shell=True, check=True)

def run_bcftools_consensus(ref:str, tmp_vcf_gz:str, outdir:str, prefix:str) -> str:
    """
    Run bcftools index
    Note: we manipulate the fasta header to conform to the format of samtools output
    """
    outfasta = os.path.join(outdir, f"{prefix}.bcftools.fasta")
    cmd = (f"bcftools consensus --fasta-ref {ref} --haplotype 1 "
           f'--missing "N" --absent "N" '
           f"{tmp_vcf_gz} "
           f"| awk '{{print $1}}' > {outfasta}"
    )
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def remove_tmp_vcf_files(tmp_vcf_gz:str, tmp_vcf_gz_index:str) -> None:
    """Remove tmp vcf files"""
    cmd = f"rm {tmp_vcf_gz} {tmp_vcf_gz_index}"
    subprocess.run(cmd, shell=True, check=True)

def remove_index_files(ref:str, bam:str) -> None:
    """Remove index files"""
    suffixes = (".amb", ".ann", ".bwt", ".fai", ".pac", ".sa")
    for suf in suffixes:
        f = ref + suf
        os.remove(f)
    bai = bam + '.bai'
    os.remove(bai)

def add_string_to_fasta_header(fasta_file:str, suf:str) -> None:
    """Add consensus method in fasta header"""
    fasta_path = Path(fasta_file)
    name = fasta_path.name
    suf = ".fasta"
    bn = name[: -len(suf)]
    ctype = bn.rsplit(".", 1)[-1]
    text = fasta_path.read_text(encoding="utf-8")
    out_lines = []
    for line in text.splitlines(True):  # keep line endings
        if line.startswith(">"):
            if line.endswith("\r\n"):
                core, nl = line[:-2], "\r\n"
            elif line.endswith("\n"):
                core, nl = line[:-1], "\n"
            else:
                core, nl = line, ""
            out_lines.append(f"{core}|{ctype}-consensus {nl}")
        else:
            out_lines.append(line)
    fasta_path.write_text("".join(out_lines), encoding="utf-8", newline="")

#def split_and_gather(ref:str, list_of_fasta:) -> None:
#    """Split reference in to separate seqs, and gather results from the different output fasta"""
    # ref is a file in fasta format containing one or several fasta headers.
    # list_of_fasta is a list of other fasta files. Each of those fasta files contain one or several fasta entries.
    # Task:
    # Open the ref fasta file. Read each fasta header and save to list named ref_headers. Read only the first string in headers:
    # For example, a header can look like this: ">OR199494.1|ITS|2  Staphylea holocarpa", and we wish to save only the part ">OR199494.1|ITS|2".
    # Iterate over each fasta file in list_of_fasta
    # For each sequence in fasta file, collect the fasta entry that matches (by grep or regex) a header in the ref_headers
    # If a match, write to a separate file named after

def main():
    """Main function."""

    args = parse_args(sys.argv[1:])

    logging.info("Start of workflow for %s", args.prefix)

    tmp_vcf_gz = os.path.join(args.outdir, '__tmp.vcf.gz')
    tmp_vcf_gz_index = os.path.join(args.outdir, '__tmp.vcf.gz.csi')
    bam = os.path.join(args.outdir, f"{args.prefix}.bam")
    list_of_fasta = []

    # Run run_bwa_index_of_ref unless file exists with suffix .amb
    ref_index = args.ref + '.amb'
    if not os.path.isfile(ref_index):
        run_bwa_index_of_ref(args.ref)

    # Run bwa mem
    logging.info("run bwa for %s using %s threads", args.prefix, args.threads)
    run_bwa(args.fastq1, args.fastq2, args.ref, args.outdir, args.prefix, args.threads)

    # Run samtools depth with -a
    logging.info("run samtools depth for %s", args.prefix)
    run_samtools_depth_a(bam, args.outdir, args.prefix, args.threads)

    # Run samtools stats
    #run_samtools_stats(bam, args.outdir, args.prefix, args.threads)

    # Run samtools index
    run_samtools_index(bam, args.threads)

    # Run samtools consensus
    logging.info("1. run samtools consensus for %s", args.prefix)
    fasta_file = run_samtools_consensus(bam, args.outdir, args.prefix, args.threads)
    add_string_to_fasta_header(fasta_file, ".fasta")
    list_of_fasta.append(fasta_file)

    # Run samtools consensus with -a option
    logging.info("2. run samtools consensus -a for %s", args.prefix)
    fasta_file = run_samtools_consensus_a(bam, args.outdir, args.prefix, args.threads)
    add_string_to_fasta_header(fasta_file, ".fasta")
    list_of_fasta.append(fasta_file)

    # Run samtools consensus with --ambig option
    logging.info("3. run samtools consensus --ambig for %s", args.prefix)
    fasta_file = run_samtools_consensus_iupac(bam, args.outdir, args.prefix, args.threads)
    add_string_to_fasta_header(fasta_file, ".fasta")
    list_of_fasta.append(fasta_file)

    # Run bcftools consensus
    logging.info("4. run bcftools consensus for %s", args.prefix)
    run_bcftools_mpileup(args.ref, args.maxdepth, bam, args.mindepth, tmp_vcf_gz)
    run_bcftools_index(tmp_vcf_gz, args.threads)
    fasta_file = run_bcftools_consensus(args.ref, tmp_vcf_gz, args.outdir, args.prefix)
    add_string_to_fasta_header(fasta_file, ".fasta")
    list_of_fasta.append(fasta_file)

    # Split reference and gather output
    #split_and_gather(args.ref, list_of_fasta)

    # Clean up temp files
    remove_index_files(args.ref, bam)
    if os.path.isfile(tmp_vcf_gz):
        remove_tmp_vcf_files(tmp_vcf_gz, tmp_vcf_gz_index)

    logging.info("End of workflow for %s", args.prefix)

if __name__ == "__main__":
    main()

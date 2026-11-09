"""Channel-wiring regression tests using lightweight stand-ins for bioinformatics tools.

Run with: python3 -m unittest discover -s tests -v
Requires Nextflow and a supported Java runtime on PATH.
"""

import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
TOOL = r"""#!/usr/bin/env python3
import gzip
from pathlib import Path
import shutil
import sys

tool = Path(sys.argv[0]).name
args = sys.argv[1:]

def value(flag):
    return args[args.index(flag) + 1]

if tool == "bwa":
    if args[0] == "index":
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            Path(args[-1] + suffix).touch()
    else:
        sys.stdout.write("@HD\tVN:1.6\n")
elif tool == "fastp":
    for read in ("1", "2"):
        with open(value("--in" + read), "rb") as source:
            with gzip.open(value("--out" + read), "wb") as target:
                shutil.copyfileobj(source, target)
    Path(value("--json")).write_text("{}")
    Path(value("--html")).write_text("<html></html>")
elif tool == "samtools":
    if args[0] == "sort":
        Path(value("-o")).write_text(sys.stdin.read())
    elif args[0] == "faidx":
        Path(args[1] + ".fai").touch()
    elif args[0] == "index":
        Path(args[-1] + ".bai").touch()
    elif args[0] == "depth":
        print("reference\t1\t1")
    elif args[0] == "consensus":
        Path(value("-o")).write_text(">reference\nACGT\n")
    else:
        raise SystemExit("Unexpected samtools arguments: " + repr(args))
elif tool == "bcftools":
    if args[0] in ("mpileup", "call"):
        if args[0] == "call":
            sys.stdin.read()
        print("mock VCF")
    elif args[0] == "view":
        Path(value("--output")).write_text(sys.stdin.read())
    elif args[0] == "index":
        Path(args[-1] + ".csi").touch()
    elif args[0] == "consensus":
        print(">reference\nACGT")
    else:
        raise SystemExit("Unexpected bcftools arguments: " + repr(args))
elif tool == "mafft":
    sys.stdout.write(Path(args[-1]).read_text())
else:
    raise SystemExit("Unexpected tool: " + tool)
"""


@unittest.skipUnless(NEXTFLOW, "Nextflow is required")
class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        tools = self.directory / "bin"
        tools.mkdir()
        executable = tools / "mock-tool"
        executable.write_text(TOOL)
        executable.chmod(0o755)
        for tool in ("bwa", "fastp", "samtools", "bcftools", "mafft"):
            (tools / tool).symlink_to(executable)
        self.env = dict(os.environ, PATH=f"{tools}{os.pathsep}{os.environ['PATH']}")
        (self.directory / "test.config").write_text(
            f"includeConfig '{ROOT / 'nextflow.config'}'\n"
            "conda.enabled = false\n"
            "process.cpus = 1\n"
            "process.maxForks = 8\n"
            "executor.queueSize = 8\n"
        )

    def fixtures(self, samples, references=1):
        refs = []
        for index in range(references):
            folder = self.directory / f"reference{index}"
            folder.mkdir()
            ref = folder / "reference.fasta"
            ref.write_text(">reference\nACGT\n")
            refs.append(ref)
        rows = []
        for index in range(samples):
            sample = f"sample{index:02}"
            for lane in range(2):
                reads = []
                for mate in (1, 2):
                    folder = self.directory / sample / f"lane{lane}"
                    folder.mkdir(parents=True, exist_ok=True)
                    read = folder / f"R{mate}.fastq.gz"
                    with gzip.open(read, "wt") as handle:
                        handle.write(f"@{sample}_lane{lane}/{mate}\nACGT\n+\nIIII\n")
                    reads.append(str(read))
                rows.append([sample, *reads, str(refs[index % references])])
        sheet = self.directory / "samples.csv"
        with sheet.open("w") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "fastq1", "fastq2", "ref"])
            writer.writerows(rows)
        return sheet, rows

    def run_workflow(self, arguments):
        result = subprocess.run(
            [NEXTFLOW, "-C", str(self.directory / "test.config"),
             "run", str(ROOT / "main.nf"), "-ansi-log", "false",
             "-with-trace", "trace.tsv", "--threads", "1",
             "--outdir", "output", *arguments],
            cwd=self.directory, env=self.env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300,
        )
        self.assertEqual(result.returncode, 0, result.stdout)
        with (self.directory / "trace.tsv").open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def check_counts(self, trace, samples, references, fastp, merge=True, mafft=False):
        expected = {
            "BWA_INDEX": references, "BWA_MEM": samples,
            "SAMTOOLS_DEPTH": samples, "SAMTOOLS_INDEX_BAM": samples,
            "SAMTOOLS_CONSENSUS": samples, "SAMTOOLS_CONSENSUS_A": samples,
            "SAMTOOLS_CONSENSUS_IUPAC": samples, "BCFTOOLS_MPILEUP_CALL": samples,
            "BCFTOOLS_INDEX": samples, "BCFTOOLS_CONSENSUS": samples,
            "CONCAT_CONSENSUS_FASTA": samples,
        }
        if merge:
            expected["MERGE_FASTQ"] = samples
        if fastp:
            expected["FASTP"] = samples
        if mafft:
            expected.update(GATHER_SEQS_FOR_MAFFT=samples, MAFFT_ALIGN=samples)
        actual = {}
        for task in trace:
            self.assertEqual(task["status"], "COMPLETED")
            process = task["name"].split(" (")[0]
            actual[process] = actual.get(process, 0) + 1
        self.assertEqual(actual, expected)
        for index in range(samples):
            fasta = self.directory / "output" / "consensus" / f"sample{index:02}.fasta"
            self.assertEqual(fasta.read_text().count(">"), 4)

    def test_55_samples_two_lanes_fastp(self):
        sheet, _ = self.fixtures(55)
        trace = self.run_workflow(["--samplesheet", str(sheet), "--fastp"])
        self.check_counts(trace, 55, 1, fastp=True)
        for index in range(55):
            sample = f"sample{index:02}"
            for mate in (1, 2):
                read = self.directory / "output" / "fastp" / f"{sample}_{mate}.trimmed.fastq.gz"
                with gzip.open(read, "rt") as handle:
                    self.assertEqual(
                        handle.read(),
                        "".join(f"@{sample}_lane{lane}/{mate}\nACGT\n+\nIIII\n" for lane in range(2)),
                    )

    def check_multiple_references(self, fastp):
        sheet, _ = self.fixtures(3, references=2)
        arguments = ["--samplesheet", str(sheet), "--mafft"]
        if fastp:
            arguments.append("--fastp")
        trace = self.run_workflow(arguments)
        self.check_counts(trace, 3, 2, fastp=fastp, mafft=True)

    def test_multiple_references_without_fastp(self):
        self.check_multiple_references(fastp=False)

    def test_multiple_references_with_fastp(self):
        self.check_multiple_references(fastp=True)

    def test_single_sample_fastp(self):
        _, rows = self.fixtures(1)
        sample, r1, r2, ref = rows[0]
        trace = self.run_workflow([
            "--prefix", sample, "--fastq1", r1, "--fastq2", r2, "--ref", ref, "--fastp",
        ])
        self.check_counts(trace, 1, 1, fastp=True, merge=False)


if __name__ == "__main__":
    unittest.main()

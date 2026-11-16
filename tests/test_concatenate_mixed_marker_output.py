"""CLI tests for selecting and concatenating marker-specific FASTA files."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "concatenate_mixed_marker_output.py"


class ConcatenateMixedMarkerOutputTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.input_folder = root / "input"
        self.output_folder = root / "output"
        self.input_folder.mkdir()
        self.filenames = (
            "plain.fasta",
            "sample.bcftools.fasta",
            "sample.samtools-a.fasta",
            "sample.samtools.fasta",
            "sample.samtools-iupac.fasta",
            "sample.samtools-a-iupac.fasta",
            "sample.custom.fa",
        )
        for filename in self.filenames:
            (self.input_folder / filename).write_text(
                f">{filename} [ref:marker|1]\nAC\nGT\n"
            )

    def run_script(self, *args):
        subprocess.run(
            [sys.executable, str(SCRIPT), str(self.input_folder),
             str(self.output_folder), *args],
            check=True,
        )
        lines = (self.output_folder / "marker_1.fasta").read_text().splitlines()
        self.assertEqual(lines[1::2], ["ACGT"] * (len(lines) // 2))
        return {line.split()[0][1:] for line in lines[::2]}

    def test_default_excludes_generated_consensus_files(self):
        self.assertEqual(self.run_script(), {"plain.fasta"})

    def test_explicit_fasta_includes_all_matching_files(self):
        self.assertEqual(self.run_script("--suffix=.fasta"), set(self.filenames[:-1]))

    def test_explicit_excluded_suffix_selects_only_matches(self):
        self.assertEqual(
            self.run_script("--suffix=.samtools-a-iupac.fasta"),
            {"sample.samtools-a-iupac.fasta"},
        )

    def test_arbitrary_suffix_matches_without_default_filter(self):
        self.assertEqual(self.run_script("--suffix=.fa"), {"sample.custom.fa"})


if __name__ == "__main__":
    unittest.main()

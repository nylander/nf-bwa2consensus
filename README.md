# nf-bwa2consensus

Map fastq using [`bwa mem`](https://github.com/lh3/BWA) against reference and
calculate consensus sequence using different methods ([`samtools
consensus`](https://github.com/samtools/samtools) and [`bcftools
consensus`](https://samtools.github.io/bcftools/howtos/consensus-sequence.html)).

Input (Illumina fastq) can optionally be filtered using
[fastp](https://github.com/opengene/fastp).

Tested using [nextflow](https://www.nextflow.io/) version 26.04.1.12112.

A python script for accomplishing basically the same task for individual
samples are located in [scripts/bwa2consensus.py](scripts/bwa2consensus.py).

## Input

1. Sample name or prefix
2. fastq files (paired-end, can be compressed)
3. Reference sequence in fasta format

These input can be specified on command line OR in a CSV file.
See [Run-examples](#run) below.

## Output

An output directory containing (Example. See also output from [single sample
run](#single-sample) below):

    output/
    ├── consensus
    │   ├── Apa.bcftools.fasta
    │   ├── Apa.fasta
    │   ├── Apa.samtools-a.fasta
    │   ├── Apa.samtools.fasta
    │   ├── Apa.samtools-iupac.fasta
    │   └── matK.fas
    └── depth
        ├── Apa.depth.tsv
        └── matK.fas

## Run

### single sample

```
$ nextflow run main.nf \
   --fastq1='input/fastq/S14_1001_01_R1.fq.gz' \
   --fastq2='input/fastq/S14_1001_01_R2.fq.gz' \
   --prefix='Apa' \
   --ref='input/reference/matK.fas' \
   --outdir='single' \
   --threads=8 \
   --mindepth=2 \
   --maxdepth=50 \
   --keepbam \
   --fastp \
   --fastp_args='--qualified_quality_phred 30 --length_required 50' \
   --mafft \
   --mafft_args='--auto'
```

Output from this single-sample run where we do fastp-filtering, keeping the bam
file, and providing a multiple-sequence alignment with reference and the
different consensus variants:

```
single/
├── bam
│   ├── Apa.bam
│   └── matK.fas
├── consensus
│   ├── Apa.bcftools.fasta
│   ├── Apa.fasta
│   ├── Apa.samtools-a.fasta
│   ├── Apa.samtools.fasta
│   ├── Apa.samtools-iupac.fasta
│   └── matK.fas
├── depth
│   ├── Apa.depth.tsv
│   └── matK.fas
├── fastp
│   ├── Apa_1.trimmed.fastq.gz
│   ├── Apa_2.trimmed.fastq.gz
│   ├── Apa.fastp.html
│   └── Apa.fastp.json
└── mafft
    └── Apa.aln.fasta
```

### three samples, different references (see [samples-matK-trnLF.csv](input/samples-matK-trnLF.csv))

```
$ nextflow run main.nf --samplesheet input/samples-matK-trnLF.csv
```

### three samples, different references in one multi-fasta file (see [samples.csv](input/samples.csv))

```
$ nextflow run main.nf --samplesheet input/samples.csv
```

For multiple read pairs per sample, repeat the sample name and reference on
each CSV row. All R1 files for a sample are merged together, and all R2 files
are merged together, before optional fastp filtering. A sample must use a
single reference file.

For example, 55 samples with two read pairs each produce 55 merged read pairs
(110 FASTQ files), 55 FASTP tasks when `--fastp` is enabled, and 55 BWA_MEM
tasks. BWA_INDEX runs once per distinct reference file, not once per sample.

**Note**: the option `--mafft` should not be used for a run with a "mixed" fasta
reference: sequences are concatenated before the MSA and can therefore create
chimeric output in this case. The output from a run "mixed fasta" run can,
however, be parsed using the script
[`concatenate_mixed_marker_output.py`](scripts/README.md).

### Remove output

    $ make clean

## Issues

The process involving `fastp` may occasionally stall the workflow. Try to abort
(`Ctrl+C`) and then run again with the `nextflow run -resume` option.

## Workflow regression tests

With Nextflow and a supported Java runtime installed, run:

```
python3 -m unittest discover -s tests -v
```

These tests run the workflow with lightweight stand-ins for the bioinformatics
tools to check channel wiring, task counts, and lane merging. They do not
validate biological results or require Conda environments.

## License

MIT [LICENSE](LICENSE)

The script [scripts/sam2consensus.py](scripts/sam2consensus.py) is under the
[GPL v.3 license](scripts/LICENSE.GPLv3).

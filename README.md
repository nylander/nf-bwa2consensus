# nf-bwa2consensus

Map fastq with [`bwa mem`](https://github.com/lh3/BWA) against reference and
calculate consensus using different methods ([`samtools
consensus`](https://github.com/samtools/samtools) and [`bcftools
consensus`](https://samtools.github.io/bcftools/howtos/consensus-sequence.html)).

Input (Illumina fastq) can optionally be filtered using
[fastp](https://github.com/opengene/fastp).

Tested using [nextflow](https://www.nextflow.io/) version 25.10.4.11173.

A python script for accomplishing basically the same task for individual samples are located in [scripts/bwa2consensus.py](scripts/bwa2consensus.py).

## Input

1. Sample name or prefix
2. fastq files (paired-end, can be compressed)
3. Reference sequence in fasta format

These input can be specified on command line OR in a CSV file.
See [Run-examples](#run) below.

## Output

An output directory containing (Example. See also output from [single sample
run](# single sample) below):

    output/
    ├── consensus
    │   ├── 28S.fas
    │   ├── S1009.bcftools.fasta
    │   ├── S1009.samtools_a.fasta
    │   ├── S1009.samtools.fasta
    │   └── S1009.samtools_iupac.fasta
    └── depth
        ├── 28S.fas
        └── S1009.depth.tsv

## Run

### single sample

```
$ nextflow run main.nf \
   --fastq1='input/fastq/P39314_1009_S17_L003_R1.fq.gz' \
   --fastq2='input/fastq/P39314_1009_S17_L003_R2.fq.gz' \
   --prefix='S1009' \
   --ref='input/reference/28S.fas' \
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

Output from this single-sample run where we do fastp-filtering, keeping the
bam file, and providing a multiple-sequence aslignment between the reference
and the different consensus variants:

```
single/
├── bam
│   ├── 28S.fas
│   └── S1009.bam
├── consensus
│   ├── 28S.fas
│   ├── S1009.bcftools.fasta
│   ├── S1009.fasta
│   ├── S1009.samtools-a.fasta
│   ├── S1009.samtools.fasta
│   └── S1009.samtools-iupac.fasta
├── depth
│   ├── 28S.fas
│   └── S1009.depth.tsv
├── fastp
│   ├── S1009_1.trimmed.fastq.gz
│   ├── S1009_2.trimmed.fastq.gz
│   ├── S1009.fastp.html
│   └── S1009.fastp.json
└── mafft
    └── S1009.aln.fasta
```

### two samples, same ref

```
$ nextflow run main.nf \
    --samplesheet input/samples.csv
```

### three samples, different refs

```
$ nextflow run main.nf \
    --samplesheet input/samples-28s-coi.csv
```

### three samples, different refs in one multi fasta

```
$ nextflow run main.nf \
    --samplesheet input/samples-mixed-ref.csv
```

### Remove output

    $ nextflow clean -f
    $ rm -rf work .nextflow* output

## License

MIT [LICENSE](LICENSE)

The script [scripts/sam2consensus.py](scripts/sam2consensus.py) is under the
[GPL v.3 license](scripts/LICENSE.GPLv3).

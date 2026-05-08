# scripts

## bwa2consensus.py

Map fastq sequences to reference using bwa mem, calculate consensus sequence
using samtools consensus and/or bcftools.

[MIT License](../LICENSE)

```
usage: bwa2consensus.py [-h] --fastq1 FASTQ1 --fastq2 FASTQ2 --ref REF [--outdir OUTDIR]
                        [--prefix PREFIX] [--threads THREADS]

Map fastq sequences to reference using bwa mem, calculate consensus sequence using samtools consensus.

options:
  -h, --help         show this help message and exit
  --fastq1 FASTQ1    Fastq file with forward reads
  --fastq2 FASTQ2    Fastq file with reverse reads
  --ref REF          Reference fasta file
  --outdir OUTDIR    Output directory
  --prefix PREFIX    Prefix for output files
  --threads THREADS  Number of threads to use

```

## sam2consensus.py

Forked from <https://github.com/edgardomortiz/sam2consensus>

[GPL v.3 license](LICENSE.GPLv3)

```
usage: sam2consensus.py [-h] -i FILENAME [-c THRESHOLDS] [-n N] [-o OUTFOLDER]
                        [-p PREFIX] [-m MIN_DEPTH] [-f FILL] [-d MAXDEL]

+------------------------------------------------------------------+
| sam2consensus.py: extract the consensus sequence from a SAM file |
+------------------------------------------------------------------+

The program takes as input a SAM file (.sam or .sam.gz) resulting from mapping
short reads to a reference (the reference sequences can correspond to separate
genes for example), then it calculates the consensus sequence from the aligned
reads alone. A single or multiple consensus thresholds can be specified, the
program also adds insertions, if many long insertions are expected, we
recommend to perform indel ralignment before for optimal results. The consensus
method is the one used by Geneious and described in detail in
http://assets.geneious.com/manual/8.1/GeneiousManualse41.html

Regions with no coverage are filled with -s (or a different character if
specified). Input SAM files don't need to be sorted. Original reference FASTAs
are not necessary since the consensus is reference-free.

It will produce a FASTA file per reference containing as many sequences as
thresholds were specified.

options:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the SAM file, SAM does not need to be sorted and
                        can be compressed with gzip
  -c THRESHOLDS, --consensus-thresholds THRESHOLDS
                        List of consensus thresold(s) separated by commas, no
                        spaces, example: -c 0.25,0.75,0.50, default=0.25
  -n N                  Split FASTA output sequences every n nucleotides,
                        default=do not split sequence
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Name of output folder, default=same folder as input
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name, default=input filename
                        without .sam extension
  -m MIN_DEPTH, --min-depth MIN_DEPTH
                        Minimum read depth at each site to report the nucleotide
                        in the consensus, default=1
  -f FILL, --fill FILL  Character for padding regions not covered in the
                        reference, default= - (gap)
  -d MAXDEL, --maxdel MAXDEL
                        Ignore deletions longer than this value, default=150
```

## concatenate_mixed_marker_output.py

Map fastq sequences to reference using bwa mem, calculate consensus sequence
using samtools consensus and/or bcftools.

Concatenate marker-specific fasta files into one fasta file per marker.
Note: fasta headers need to have labels following a specific format. For
example: `>... [ref:PV424158.1|ITS|1]`.

[MIT License](../LICENSE)

```
usage: concatenate_mixed_marker_output.py [-h] input_folder output_folder

Concatenate marker-specific fasta files into one fasta file per marker.
Note: fasta headers need to have labels following a specific format. For
example: >... [ref:PV424158.1|ITS|1]

positional arguments:
  input_folder   Path to the input folder containing fasta files
  output_folder  Path to the output folder where fasta files will be saved

options:
  -h, --help     show this help message and exit
```

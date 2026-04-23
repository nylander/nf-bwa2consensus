#
# Makefile for nf-bwa2consensus
#

.PHONY: all run clean distclean test mafft fastp fastp-mafft-keepbam

all: run

run:
	nextflow run main.nf --samplesheet input/samples.csv

fastp:
	nextflow run main.nf --samplesheet input/samples.csv --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50'

mafft:
	nextflow run main.nf --samplesheet input/samples.csv --mafft --mafft_args='--auto'

fastp-mafft-keepbam:
	nextflow run main.nf --samplesheet input/samples.csv --mafft --mafft_args='--auto' --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50' --keepbam

test: fastp-mafft-keepbam

clean:
	rm -rf work/ .nextflow*

distclean:
	rm -rf work/ output/ .nextflow*

# vim:ft=make
#

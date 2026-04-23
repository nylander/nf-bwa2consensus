#
# Makefile for nf-bwa2consensus
#

.PHONY: all run clean distclean test mafft fastp fastp-mafft-keepbam

all: run

run:
	@echo "map S1009 and S1011 against 28S"
	nextflow run main.nf --samplesheet input/samples.csv
	@tree output

two-three:
	@echo "S1009 mapped against 28S, S1010 and S1011 against COI"
	nextflow run main.nf --samplesheet input/samples-28s-coi.csv
	@tree output

fastp:
	@echo "run fastp and then map S1009 and S1011 against 28S"
	nextflow run main.nf --samplesheet input/samples.csv --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50'
	@tree output

mafft:
	@echo "map S1009 and S1011 against 28S, align ref+consensus with mafft"
	nextflow run main.nf --samplesheet input/samples.csv --mafft --mafft_args='--auto'
	@tree output

fastp-mafft-keepbam:
	@echo "run all options"
	nextflow run main.nf --samplesheet input/samples.csv --mafft --mafft_args='--auto' --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50' --keepbam
	@tree output

test: fastp-mafft-keepbam

clean:
	rm -rf work/ .nextflow*

distclean:
	rm -rf work/ output/ .nextflow*

# vim:ft=make
#

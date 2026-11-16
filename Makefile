#
# Makefile for nf-bwa2consensus
#

.PHONY: all run clean distclean test mafft keepbam fastp fastp-mafft-keepbam

all: run

run:
	@echo "map all against mix"
	nextflow run main.nf --samplesheet input/samples.csv
	@tree output

two-genes:
	@echo "map Apa against matK, Bpa,Cpa,Dpa against trnLF"
	nextflow run main.nf --samplesheet input/samples-matK-trnLF.csv
	@tree output

keepbam:
	@echo "map Apa against matK, Bpa,Cpa,Dpa against trnLF, keep bam files"
	nextflow run main.nf --samplesheet input/samples-matK-trnLF.csv --keepbam
	@tree output

fastp:
	@echo "run fastp and then map all against mix"
	nextflow run main.nf --samplesheet input/samples.csv --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50'
	@tree output

mafft:
	@echo "map Apa against matK, Bpa,Cpa,Dpa against trnLF, align ref+consensus with mafft"
	nextflow run main.nf --samplesheet input/samples-matK-trnLF.csv --mafft --mafft_args='--auto'
	@tree output

fastp-mafft-keepbam:
	@echo "run all options"
	nextflow run main.nf --samplesheet input/samples-matK-trnLF.csv --mafft --mafft_args='--auto' --fastp --fastp_args='--qualified_quality_phred 30 --length_required 50' --bwa_args='-k 15' --bwa_index_args='-a is' --keepbam
	@tree output

test: distclean fastp-mafft-keepbam

regression-test:
	python3 -m unittest discover -s tests -v

clean:
	rm -rf work/ .nextflow*

distclean:
	rm -rf work/ output/ .nextflow* tests/__pycache__/

# vim:ft=make
#

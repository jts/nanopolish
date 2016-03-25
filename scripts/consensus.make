SHELL=/bin/bash -o pipefail
.SECONDARY:

BWA=bwa

# Get the path to the Makefile
ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# Get the name of the input file without a suffix
READS_BASE=$(basename $(READS))

#
# A pipeline to recompute a consensus sequence for an assembly
#
all: $(READS_BASE).pp.sorted.bam $(READS_BASE).pp.sorted.bam.bai $(ASSEMBLY).fai

#
# Preprocess the reads to make a name map
# and uniquify names
#
%.pp.fa: %.fa
	$(ROOT_DIR)/consensus-preprocess.pl $^ > $@

# handle .fasta too
%.pp.fa: %.fasta
	$(ROOT_DIR)/consensus-preprocess.pl $^ > $@

#
# Make bwa index files for the assembly
#
$(ASSEMBLY).bwt: $(ASSEMBLY)
	$(BWA) index $^

#
# Map reads to the assembly using bwa
#
%.pp.bam: %.pp.fa $(ASSEMBLY).bwt
	$(BWA) mem -x ont2d -t 8 $(ASSEMBLY) $< | samtools view -Sb - > $@

#
# Sort BAM
#
%.sorted.bam: %.bam
	samtools sort -f $^ $@

#
# Index BAM
#
%.sorted.bam.bai: %.sorted.bam
	samtools index $^

#
# Index assembly
#
$(ASSEMBLY).fai: $(ASSEMBLY)
	samtools faidx $<

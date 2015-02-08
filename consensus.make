SHELL=/bin/bash -o pipefail
.SECONDARY:

BWA=bwa

#
# A pipeline to recompute a consensus sequence for an assembly
#
all=$(READS).pp.sorted.bam $(READS).pp.sorted.bam.bai

#
# Preprocess the reads to make a name map
# and uniquify names
#
%.pp.fa: %.fa
	./consensus-preprocess.pl $^ > $@

#
# Make bwa index files for the assembly
#
$(ASSEMBLY).bwt: $(ASSEMBLY)
	$(BWA) index $^

#
# Map reads to the assembly using bwa
#
%.pp.bam: %.pp.fa $(ASSEMBLY).bwt
	$(BWA) mem -x ont2d -t 8 $(ASSEMBLY) $< > $@

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

SHELL=/bin/bash -o pipefail
.PRECIOUS: $(ASSEMBLY).ranges $(ASSEMBLY).bwt
.SECONDARY:

BWA=bwa
THREADS=8

TARGETS = $(ASSEMBLY)-$(READS).pp.sorted.bam $(ASSEMBLY)-$(READS).pp.sorted.bam.bai $(READS).pp.fa $(READS).pp.fa.fast5.fofn \
	nanopolish-$(ASSEMBLY)-$(READS).results nanopolish-$(ASSEMBLY)-$(READS).polished.fa


#
# A pipeline to recompute a consensus sequence for an assembly
#
all: $(TARGETS)

#
# Preprocess the reads to make a name map
# and uniquify names
#
%.pp.fa: %
	./consensus-preprocess.pl $^ > $@.tmp && mv $@.tmp $@

%.pp.fa.fast5.fofn: % %.pp.fa
	mv $<.fast5.fofn $@

#
# Make bwa index files for the assembly
#
$(ASSEMBLY).bwt: $(ASSEMBLY)
	$(BWA) index $^

#
# Map reads to the assembly using bwa
#
$(ASSEMBLY)-%.pp.bam: %.pp.fa $(ASSEMBLY).bwt
	$(BWA) mem -x ont2d -t $(THREADS) $(ASSEMBLY) $< | samtools view -Sb - > $@.tmp && mv $@.tmp $@

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

%.ranges : %
	python nanopolish_makerange.py $^ > $@.tmp && mv $@.tmp $@

nanopolish-$(ASSEMBLY)-$(READS).results : $(ASSEMBLY).ranges $(ASSEMBLY)-$(READS).pp.sorted.bam $(ASSEMBLY)-$(READS).pp.sorted.bam.bai $(ASSEMBLY) $(READS).pp.fa $(READS).pp.fa.fast5.fofn
	parallel --results $@.tmp -P $(THREADS) ./nanopolish consensus -o nanopolish-$(ASSEMBLY)-$(READS).{1}.fa -w {1} --r $(READS).pp.fa -b $(ASSEMBLY)-$(READS).pp.sorted.bam -g $(ASSEMBLY) -t 4 < $(ASSEMBLY).ranges \
	&& mv $@.tmp $@

nanopolish-$(ASSEMBLY)-$(READS).polished.fa  : nanopolish-$(ASSEMBLY)-$(READS).results
	python nanopolish_merge.py $(ASSEMBLY) nanopolish-$(ASSEMBLY)-$(READS).*.fa > $@.tmp && mv $@.tmp $@

clean:
	rm -rf $(TARGETS) *.tmp

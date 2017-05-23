# Nanopolish

[![Build Status](https://travis-ci.org/jts/nanopolish.svg?branch=master)](https://travis-ci.org/jts/nanopolish)

Software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish is designed to calculate a new consensus sequence for a draft genome assembly produced by a program like [CANU](https://github.com/marbl/canu). It does not error correct reads.

## Dependencies

[libhdf5](http://www.hdfgroup.org/HDF5/release/obtain5.html). It is automatically downloaded and compiled during `make` step but you can disable it with: `HDF5=nofetch make`. It is not necessary to install it (and `make install` is not called either). The nanopolish binary will link using a libhdf5.a (statically).

[eigen](http://eigen.tuxfamily.org). It is automatically downloaded and compiled/ Currently you cannot override that.

[biopython](http://www.biopython.org)

A compiler that supports C++11 is need to build the sources. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/). libhdf5 can be automatically installed by the Makefile if you do not have it already (see below).

## Installation instructions

You will need to run ```git clone --recursive https://github.com/jts/nanopolish.git``` to get the source code and submodules. You can then compile nanopolish by running:

```
make
```

This will automatically download and install libhdf5.

## Nanopolish modules

The main subprograms of nanopolish are:

```
nanopolish extract: extract reads in FASTA or FASTQ format from a directory of FAST5 files
nanopolish call-methylation: predict genomic bases that may be methylated
nanopolish variants: detect SNPs and indels with respect to a reference genome
nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly
nanopolish eventalign: align signal-level events to k-mers of a reference genome
```

## Analysis workflow examples

### Computing a new consensus sequence for a draft assembly

First we prepare the data by extracting the reads from the FAST5 files, and aligning the basecalls to our draft assembly (`draft.fa`).

```
# Extract the QC-passed reads from a directory of FAST5 files
nanopolish extract --type [2d|template] directory/pass/ > reads.fa

# Index the draft genome (produced by [CANU](https://github.com/marbl/canu) with megabase-sized contigs)
bwa index draft.fa

# Align the basecalled reads to the draft sequence
bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam
```

Now, we use nanopolish to compute the consensus sequence (50kb blocks of the genome will be output). We'll run this in parallel:

```
python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1
```

This command will run the consensus algorithm on eight 50kbp segments of the genome at a time, using 4 threads each. Change the ```-P``` and ```--threads``` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual 50kb segments together back into the final assembly:

```
python nanopolish_merge.py polished.*.fa > polished_genome.fa
```

## Calling Methylation

nanopolish can use the signal-level information measured by the sequencer to detect 5-mC as described [here](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4184.html). Here's how you run it:

```
# Extract all reads from a directory of FAST5 files
nanopolish extract -r --type template directory/ > reads.fa

# Align the basecalled reads to a reference genome
bwa mem -x ont2d -t 8 reference.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam

# Run the methylation caller
nanopolish call-methylation -t 8 -r reads.fa -g reference.fa -b reads.sorted.bam > methylation.tsv
```

The output of call-methylation is a tab-separated file containing per-read log-likelihood ratios (positive values indicate more evidence for 5-mC, negative values indicate more evidence for C). Each line contains the name of the read that covered the CpG site, which allows methylation calls to be phased. We have provided a script to calculate per-site methylation frequencies using call-methylation's output:

```
python /path/to/nanopolish/scripts/calculate_methylation_frequency -c 2.5 -i methylation.tsv > frequencies.tsv
```

The output of this script is a tab-seperated file containing the genomic position of the CpG site, the number of reads that covered the site, and the percentage of those reads that were predicted to be methylated. The `-c 2.5` option requires the absolute value of the log-likelihood ratio to be at least 2.5 to make a call, otherwise the read will be ignored. This helps reduce calling errors as only sites with sufficient evidence will be included in the calculation.

## To run using docker

First build the image from the dockerfile:
```
docker build .
```
Note the uuid given upon successful build.
Then you can run nanopolish from the image:
```
docker run -v /path/to/local/data/data/:/data/ -it :image_id  ./nanopolish eventalign -r /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
```

## Credits and Thanks

The fast table-driven logsum implementation was provided by Sean Eddy as public domain code. This code was originally part of [hmmer3](http://hmmer.janelia.org/).

# Nanopolish

[![Build Status](https://travis-ci.org/jts/nanopolish.svg?branch=master)](https://travis-ci.org/jts/nanopolish)

Software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish can calculate an improved consensus sequence for a draft genome assembly, detect base modifications, call SNPs and indels with respect to a reference genome and more (see Nanopolish modules, below).

## Dependencies

[libhdf5](http://www.hdfgroup.org/HDF5/release/obtain5.html) is automatically downloaded and compiled when running `make` but this can be disabled with: `HDF5=nofetch make`. The nanopolish binary will link libhdf5.a statically.

[eigen](http://eigen.tuxfamily.org) is also automatically downloaded and included when compiling with `make`.

[biopython](http://www.biopython.org) is required to run the helpers in `scripts/`.

[htslib](http://github.com/samtools/htslib) is included as a submodule and compiled automatically.

A compiler that supports C++11 is needed to build the sources. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/).

## Installation instructions

### Installing the latest code from github (recommended)

You can download and compile the latest code from github as follows:

```
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
make
```

### Installing a particular release

When major features have been added or bugs fixed, we will tag and release a new version of nanopolish. If you wish to use a particular version, you can checkout the tagged version before compiling:

```
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
git checkout v0.7.1
make
```

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

### Data preprocessing

Nanopolish needs access to the signal-level data measured by the nanopore sequencer. The first step of any nanopolish workflow is to prepare the input data by telling nanopolish where to find the signal files. If you ran Albacore 2.0 on your data you should run `nanopolish index` on your input reads (-d can be specified more than once if using multiple runs):

```
# Only run this if you used Albacore 2.0 or later
nanopolish index -d /path/to/raw_fast5s/ albacore_output.fastq
```

If you basecalled your reads with Albacore 1.2 or earlier, you should run `nanopolish extract` on your input reads instead:

```
# Only run this if you used Albacore 1.2 or later
nanopolish extract --type template directory/pass/ > reads.fa
```

Note these two commands are mutually exclusive - you only need to run one of them. You need to decide what command to run depending on the version of Albacore that you used. In the following sections we assume you have preprocessed the data by following the instructions above and that your reads are in a file named `reads.fa`.

### Computing a new consensus sequence for a draft assembly

The original purpose of nanopolish was to compute an improved consensus sequence for a draft genome assembly produced by a long-read assembly like [canu](https://github.com/marbl/canu). This section describes how to do this, starting with your draft assembly which should have megabase-sized contigs.

```
# Index the draft genome
bwa index draft.fa

# Align the basecalled reads to the draft sequence
bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam
```

Now, we use nanopolish to compute the consensus sequence (the genome is polished in 50kb blocks and there will be one output file per block). We'll run this in parallel:

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

The fast table-driven logsum implementation was provided by Sean Eddy as public domain code. This code was originally part of [hmmer3](http://hmmer.janelia.org/). Nanopolish also includes code from Oxford Nanopore's [scrappie](https://github.com/nanoporetech/scrappie) basecaller. This code is licensed under the MPL.

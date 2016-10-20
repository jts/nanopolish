# Nanopolish

[![Build Status](https://travis-ci.org/jts/nanopolish.svg?branch=master)](https://travis-ci.org/jts/nanopolish)

A nanopore consensus algorithm using a signal-level hidden Markov model.

## Dependencies

The program requires [libhdf5](http://www.hdfgroup.org/HDF5/release/obtain5.html) and a compiler that supports C++11. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/). libhdf5 can be automatically installed by the Makefile if you do not have it already (see below).

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
nanopolish eventalign: align signal-level events to k-mers of a reference genome
nanopolish variants: detect SNPs and indels with respect to a reference genome
nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly
```

## Analysis workflows

The two main uses of nanopolish are to calculate an improved consensus sequence for a draft genome assembly, and to find SNPs and indels with respect to a reference genome.

### Computing a new consensus sequence for a draft assembly

First we prepare the data by extracting the reads from the FAST5 files, and aligning them in base and event space to our draft assembly (`draft.fa`).

```
# Extract the QC-passed reads from a directory of FAST5 files
nanopolish extract --type [2d|template] directory/pass/ > reads.fa

# Index the reference genome
bwa index draft.fa

# Align the reads in base space
bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools view -Sb - | samtools sort -f - reads.sorted.bam
samtools index reads.sorted.bam

# Copy the nanopolish model files into the working directory
cp /path/to/nanopolish/etc/r9-models/* .

# Align the reads in event space
nanopolish eventalign -t 8 --sam -r reads.fa -b reads.sorted.bam -g draft.fa --models nanopolish_models.fofn | samtools view -Sb - | samtools sort -f - reads.eventalign.sorted.bam
samtools index reads.eventalign.sorted.bam
```

Now, we use nanopolish to compute the consensus sequence. We'll run this in parallel:

```
python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -e reads.eventalign.sorted.bam -t 4 --min-candidate-frequency 0.1 --models nanopolish_models.fofn
```

This command will run the consensus algorithm on eight 10kbp segments of the genome at a time, using 4 threads each. Change the ```-P``` and ```--threads``` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual segments together into the final assembly:

```
python nanopolish_merge.py polished.*.fa > polished_genome.fa
```

## Fixing homopolymers

Nanopolish 0.5 contains an experimental ```--fix-homopolymers``` option that will use event durations to improve the consensus accuracy around homopolymers. This option has only been tested on deep (>100X) data where it gives a minor improvement in accuracy. It is left off by default for now until it is tested further.

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

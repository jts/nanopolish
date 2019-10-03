# Nanopolish

[![Build Status](https://travis-ci.org/jts/nanopolish.svg?branch=master)](https://travis-ci.org/jts/nanopolish)

Software package for signal-level analysis of Oxford Nanopore sequencing data. Nanopolish can calculate an improved consensus sequence for a draft genome assembly, detect base modifications, call SNPs and indels with respect to a reference genome and more (see Nanopolish modules, below).

## Release notes

* 0.11.1: `nanopolish polya` now supports SQK-RNA-002 kits with automatic backwards-compatibility with SQK-RNA-001

* 0.11.0: support for multi-fast5 files. `nanopolish methyltrain` now subsamples input data, improving speed and memory usage

* 0.10.2: added new program `nanopolish polya` to estimate the length of poly-A tails on direct RNA reads (by @paultsw)

* 0.10.1: `nanopolish variants --consensus` now only outputs a VCF file instead of a fasta sequence. The VCF file describes the changes that need to be made to turn the draft sequence into the polished assembly. A new program, `nanopolish vcf2fasta`, is provided to generate the polished genome (this replaces `nanopolish_merge.py`, see usage instructions below). This change is to avoid issues when merging segments that end on repeat boundaries (reported by Michael Wykes and Chris Wright).

## Dependencies

A compiler that supports C++11 is needed to build nanopolish. Development of the code is performed using [gcc-4.8](https://gcc.gnu.org/gcc-4.8/).

By default, nanopolish will download and compile all of its required dependencies. Some users however may want to use system-wide versions of the libraries. To turn off the automatic installation of dependencies set `HDF5=noinstall`, `EIGEN=noinstall` or `HTS=noinstall` parameters when running `make` as appropriate. The current versions and compile options for the dependencies are:

* [libhdf5-1.8.14](http://www.hdfgroup.org/HDF5/release/obtain5.html) compiled with multi-threading support `--enable-threadsafe`
* [eigen-3.2.5](http://eigen.tuxfamily.org)
* [htslib-1.4](http://github.com/samtools/htslib) 

Additionally the helper `scripts` require [biopython](http://www.biopython.org) and [pysam](http://pysam.readthedocs.io/en/latest/installation.html).


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
git checkout v0.9.2
make
```

## Nanopolish modules

The main subprograms of nanopolish are:

```
nanopolish call-methylation: predict genomic bases that may be methylated
nanopolish variants: detect SNPs and indels with respect to a reference genome
nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly
nanopolish eventalign: align signal-level events to k-mers of a reference genome
```

## Analysis workflow examples

### Data preprocessing

Nanopolish needs access to the signal-level data measured by the nanopore sequencer. The first step of any nanopolish workflow is to prepare the input data by telling nanopolish where to find the signal files. If you ran Albacore 2.0 on your data you should run `nanopolish index` on your input reads (-d can be specified more than once if using multiple runs):

```
# Index the output of the albacore basecaller
nanopolish index -d /path/to/raw_fast5s/ -s sequencing_summary.txt albacore_output.fastq
```

The `-s` option tells nanopolish to read the `sequencing_summary.txt` file from Albacore to speed up indexing. Without this option `nanopolish index` is extremely slow as it needs to read every fast5 file individually. If you basecalled your run in parallel, so you have multiple `sequencing_summary.txt` files, you can use the `-f` option to pass in a file containing the paths to the sequencing summary files (one per line).

### Computing a new consensus sequence for a draft assembly

The original purpose of nanopolish was to compute an improved consensus sequence for a draft genome assembly produced by a long-read assembly like [canu](https://github.com/marbl/canu). This section describes how to do this, starting with your draft assembly which should have megabase-sized contigs. We've also posted a tutorial including example data [here](http://nanopolish.readthedocs.io/en/latest/quickstart_consensus.html).

```
# Index the draft genome
bwa index draft.fa

# Align the basecalled reads to the draft sequence
bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
samtools index reads.sorted.bam
```

Now, we use nanopolish to compute the consensus sequence (the genome is polished in 50kb blocks and there will be one output file per block). We'll run this in parallel:

```
python3 nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1
```

This command will run the consensus algorithm on eight 50kbp segments of the genome at a time, using 4 threads each. Change the ```-P``` and ```--threads``` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual 50kb segments together back into the final assembly:

```
nanopolish vcf2fasta -g draft.fa polished.*.vcf > polished_genome.fa
```

## Calling Methylation

nanopolish can use the signal-level information measured by the sequencer to detect 5-mC as described [here](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4184.html). We've posted a tutorial on how to call methylation [here](http://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html).

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

## GPU acceleration

The nanopolish consensus improvement algorithm can be performed faster using CUDA-enabled GPU acceleration. This is an experimental feature, to try this feature run with the `--gpu=1` flag e.g:
```
nanopolish variants --consensus polished_gpu.fa -w "tig00000001:200000-230000" -r reads.fasta -b reads.sorted.bam -g draft.fa --threads=8 --gpu=1
```
Note that this feature requires nanopolish to be compiled with `make cuda=1`. You should have the [CUDA toolkit installed and configured](https://docs.nvidia.com/cuda/cuda-quick-start-guide/). If your CUDA installation is not in the default location, you can provide the path to make as `make cuda=1 NVCC=/path/to/nvidia_c_compiler CUDA_LIB=/path/to/cuda/lib CUDA_INCLUDE=/path/to/cuda/include`.

## Credits and Thanks

The fast table-driven logsum implementation was provided by Sean Eddy as public domain code. This code was originally part of [hmmer3](http://hmmer.janelia.org/). Nanopolish also includes code from Oxford Nanopore's [scrappie](https://github.com/nanoporetech/scrappie) basecaller. This code is licensed under the MPL. GPU support was originally developed by Mike Vella. Hasindu Gamaarachchi revised Mike's implementation and integrated it into the main nanopolish branch. 

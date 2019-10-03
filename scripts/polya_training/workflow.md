Retraining the poly(A) segmentation hidden markov model
=======================================================
_Author: Paul Tang (paul.tang@oicr.on.ca)_

This document provides a guide to the (semi-automated) process of retraining the hidden markov model at
the core of the `nanopolish polya` subprogram.

I. Initial setup
----------------
The initial preprocessing stages for the files are the same as the typical workflow for calling poly(A) tail lengths.
Assuming you are starting off with multi-read FAST5s and basecalled FASTQs from Guppy, perform the following steps:

1. Split multi-read FAST5 files into single-read FAST5 files using `ont-fast5-api`.
2. `cat` all the FASTQ files together into a single FASTQ.
3. Generate SAM files with `minimap2` in splice-aware mode:
  `minimap2 -ax splice -uf -k14 <REF>.fa <READS>.fq > <ALNS>.sam`
4. Convert to BAM, sort, and index with `samtools {sort, index}`.

II. Iterative model development loop
------------------------------------
It is suggested to have a small development set of well-chosen reads not larger than ~500 examples; you want this to
be a large enough dataset where aggregate statistics are meaningful while also being able to quickly re-run the poly(A)
length estimation process.

As a practical note, I like to have the following `tmux` configuration across 4 panes:
* pane 1: a small pane for repeatedly running `nanopolish polya` on the dataset in question.
* pane 2: a small pane for recompiling nanopolish (i.e., re-running `make` in the top-level `nanopolish` directory).
* pane 3: a large pane (maybe 2/3 of the left or right half of your terminal screen) for making edits to the nanopolish
  source file that implements the segmentation HMM (`nanopolish/src/nanopolish_polya_estimator.cpp`).
* pane 4: a large pane (perhaps the left or right half of your terminal screen) for running the python-based interactive
  visualization/dumping/fitting scripts.

Roughly, my suggested model iteration workflow consists of the following:
1. Run `nanopolish polya --reads={...} --bam={...} --genome={...} > {...}.polya.tsv` to generate new segmentations.
2. Visualize the segmentations using the python script (see `hmmplot`).
3. Look for errors in the segmentation and modify the `nanopolish_polya_estimator.cpp` file accordingly; recompile.
4. `GOTO 1`.

In step (3), it may be helpful to obtain the raw picoamp signals for each segment and then retrain gaussian mixture
models to fit them. To do this, we recommend the following steps:
1. select a number of read ids representing segmentations visually observed to be "correct" and re-run the `polya`
   workflow on this subset;
2. generate verbose poly(A) calls with `nanopolish polya -v {...}`;
3. use `dump_signal.py` to generate an HDF5 in the desired directory containing normalized pico-ampere samples;
4. use `retrain_gmm.py` to estimate the parameters of a gaussian mixture model for each read.

III. Description of included scripts
------------------------------------
This directory contains the core executable plotting script for inspecting segmentations, as well as a script
for dumping scaled sample values and a script for easy retraining of mixture models.

The following files are included:
* `hmmplot.py`: a script for visualizing segmentations, overlaid as color-coded regions on a signal trace
  of the underlying pico-Ampere current of the read.
* `dump_signal.py`: takes a poly(A) output TSV generated with the `-v` verbose flag and dumps an HDF5 file
  containing the normalized picoampere samples for each segment of the reads of interest.
* `retrain_emission.py`: take the arrays from a run of `dump_signal.py` and output the coefficients and parameters
  of a gaussian mixture model for each segment.

Run `python <SCRIPT> --help` for further usage instructions for each script.

The files above required a suite of numerical packages to be installed as dependencies. A full list of dependencies,
in a conda-compatible environment file, is available at `hmmplot/environment.yml`.
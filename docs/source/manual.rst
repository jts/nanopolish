.. _manual:

Manual
===================

Modules available: ::

    nanopolish extract: extract reads in FASTA or FASTQ format from a directory of FAST5 files
    nanopolish call-methylation: predict genomic bases that may be methylated
    nanopolish variants: detect SNPs and indels with respect to a reference genome
    nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly
    nanopolish eventalign: align signal-level events to k-mers of a reference genome

|

Extract
--------------------

Overview
"""""""""""""""""""""""

This module is used to extract reads in FASTA or FASTQ format from a directory of FAST5 files.  

Input
"""""""""""""""""""""""

    * path to a directory of FAST5 files modified to contain basecall information

Output
"""""""""""""""""""""""

    * sequences of reads in FASTA or FASTQ format

Usage example
"""""""""""""""""""""""

::

   nanopolish extract [OPTIONS] <fast5|dir>

.. list-table:: Command line arguments
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * -  <fast5|dir>
     - Y
     - NA
     - FAST5 or path to directory of FAST5 files.

   * - -r, --recurse
     - N
     - NA
     - Recurse into subdirectories

   * - -q, --fastq
     - N
     - fasta format
     - Use when you want to extract to FASTQ format

   * - -t, --type=TYPE
     - N
     - 2d-or-template
     - The type of read either: {template, complement, 2d, 2d-or-template, any}

   * - -b, --basecaller=NAME[:VERSION]
     - N
     - NA
     - consider only data produced by basecaller NAME, optionally with given exact VERSION

   * - -o, --output=FILE
     - N
     - stdout
     - Write output to FILE

index
--------------------

Overview
"""""""""""""""""""""""
Build an index mapping from basecalled reads to the signals measured by the sequencer

Input
""""""""
    * path to directory of raw nanopore sequencing data in FAST5 format
    * basecalled reads

Output
""""""""
    * gzipped FASTA format of basecalled reads
    * index files (fai, gzi, readdb)

Readdb file format
""""""""""""""""""""
Readdb file is a tab-separated file that contains two columns. One column represents read ids and the other column represents the corresponding path to FA
ST5 file: ::

    read_id_1   /path/to/fast5/containing/reads_id_1/signals
    read_id_2   /path/to/fast5/containing/read_id_2/signals

Usage example
""""""""""""""
::

    nanopolish index [OPTIONS] -d nanopore_raw_file_directory reads.fastq

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - -d, --directory
     - Y
     - NA
     - FAST5 or path to directory of FAST5 files containing ONT sequencing raw signal information.

   * - -f, --fast5-fofn
     - N
     - NA
     - file containing the paths to each fast5 for the run



call-methylation
--------------------

Overview
"""""""""""""""""""""""

Classify nucleotides as methylated or not.

Input
"""""""""""""""""""""""

    * Basecalled ONT reads in FASTA format

Output
"""""""""""""""""""""""

    * tab-separated file containing per-read log-likelihood ratios

Usage example
"""""""""""""""""""""""

::

   nanopolish call-[OPTIONS] <fast5|dir>

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - -r, --reads=FILE
     - Y
     - NA
     - the 2D ONT reads are in fasta FILE

   * - -b, --bam=FILE
     - Y
     - NA 
     - the reads aligned to the genome assembly are in bam FILE

   * - -g, --genome=FILE
     - Y
     - NA 
     - the genome we are computing a consensus for is in FILE

   * - -t, --threads=NUM
     - N
     - 1
     - use NUM threads

   * - --progress
     - N
     - NA
     - print out a progress message

variants
--------------------

Overview
"""""""""""""""""""""""

This module is used to call single nucleotide polymorphisms (SNPs) using a signal-level HMM.  

Input
"""""""""""""""""""""""

    * basecalled reads
    * alignment info
    * genome assembly

Output
"""""""""""""""""""

    * 

Usage example
"""""""""""""""""""""""

::

   nanopolish extract [OPTIONS] <fast5|dir>

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - --snps
   * - N
   * - NA
   * - use flag to only call SNPs

   * - --consensus=FILE
   * - N
   * - NA
   * - run in consensus calling mode and write polished sequence to FILE

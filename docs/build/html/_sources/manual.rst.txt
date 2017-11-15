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

extract
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

.. list-table:: Command-line arguments
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

call-methylation
--------------------


variants
--------------------

variants --consensus
--------------------

eventalign
--------------------

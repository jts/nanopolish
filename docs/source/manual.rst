.. _manual:

Manual
===================

Modules available: ::

    nanopolish extract: extract reads in FASTA or FASTQ format from a directory of FAST5 files
    nanopolish call-methylation: predict genomic bases that may be methylated
    nanopolish variants: detect SNPs and indels with respect to a reference genome
    nanopolish variants --consensus: calculate an improved consensus sequence for a draft genome assembly
    nanopolish eventalign: align signal-level events to k-mers of a reference genome
    nanopolish phase-reads: Phase reads using heterozygous SNVs with respect to a reference genome 
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

.. list-table:: 
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

   * - ``-r``, ``--recurse``
     - N
     - NA
     - Recurse into subdirectories

   * - ``-q``, ``--fastq``
     - N
     - fasta format
     - Use when you want to extract to FASTQ format

   * - ``-t``, ``--type=TYPE``
     - N
     - 2d-or-template
     - The type of read either: {template, complement, 2d, 2d-or-template, any}

   * - ``-b``, ``--basecaller=NAME[:VERSION]``
     - N
     - NA
     - consider only data produced by basecaller NAME, optionally with given exact VERSION

   * - ``-o``, ``--output=FILE``
     - N
     - stdout
     - Write output to FILE

|

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
Readdb file is a tab-separated file that contains two columns. One column represents read ids and the other column represents the corresponding path to FAST5 file: ::

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

   * - ``-d``, ``--directory``
     - Y
     - NA
     - FAST5 or path to directory of FAST5 files containing ONT sequencing raw signal information.

   * - ``-f``, ``--fast5-fofn``
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

   nanopolish call-methylation [OPTIONS] <fast5|dir>

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``-r``, ``--reads=FILE``
     - Y
     - NA
     - the ONT reads are in fasta FILE

   * - ``-b``, ``--bam=FILE``
     - Y
     - NA 
     - the reads aligned to the genome assembly are in bam FILE

   * - ``-g``, ``--genome=FILE``
     - Y
     - NA 
     - the genome we are computing a consensus for is in FILE

   * - ``-t``, ``--threads=NUM``
     - N
     - 1
     - use NUM threads

   * - ``--progress``
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

    * VCF file

Usage example
"""""""""""""""""""""""

::

   nanopolish variants [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``--snps``
     - N
     - NA
     - use flag to only call SNPs

   * - ``--consensus=FILE``
     - N
     - NA
     - run in consensus calling mode and write polished sequence to FILE

   * - ``--fix-homopolymers``
     - N
     - NA
     - use flag to run the experimental homopolymer caller

   * - ``--faster``
     - N
     - NA
     - minimize compute time while slightly reducing consensus accuracy

   * - ``-w``, ``--window=STR``
     - N
     - NA
     - find variants in window STR (format: <chromsome_name>:<start>-<end>)

   * - ``-r``, ``--reads=FILE``
     - Y
     - NA
     - the ONT reads are in fasta FILE

   * - ``-b``, ``--bam=FILE``
     - Y
     - NA
     - the reads aligned to the reference genome are in bam FILE 

   * - ``-e``, ``--event-bam=FILE``
     - Y
     - NA
     - the events aligned to the reference genome are in bam FILE

   * - ``-g``, ``--genome=FILE``
     - Y
     - NA
     - the reference genome is in FILE

   * - ``-o``, ``--outfile=FILE``
     - N
     - stdout
     - write result to FILE

   * - ``-t``, ``--threads=NUM``
     - N
     - 1
     - use NUM threads

   * - ``-m``, ``--min-candidate-frequency=F``
     - N
     - 0.2
     - extract candidate variants from the aligned reads when the variant frequency is at least F

   * - ``-d``, ``--min-candidate-depth=D``
     - N
     - 20
     - extract candidate variants from the aligned reads when the depth is at least D

   * - ``-x``, ``--max-haplotypes=N``
     - N
     - 1000
     - consider at most N haplotypes combinations

   * - ``--max-rounds=N``
     - N
     - 50
     - perform N rounds of consensus sequence improvement

   * - ``-c``, ``--candidates=VCF``
     - N
     - NA
     - read variants candidates from VCF, rather than discovering them from aligned reads

   * - ``-a``, ``--alternative-basecalls-bam=FILE``
     - N
     - NA
     - if an alternative basecaller was used that does not output event annotations then use basecalled sequences from FILE. The signal-level events will still be taken from the -b bam

   * - ``--calculate-all-support``
     - N
     - NA
     - when making a call, also calculate the support of the 3 other possible bases

   * - ``--models-fofn=FILE``
     - N
     - NA
     - read alternatives k-mer models from FILE


event align
--------------------

Overview
"""""""""""""""""""""""

Align nanopore events to reference k-mers

Input
"""""""""""""""""""""""

    * basecalled reads
    * alignment information
    * assembled genome

Usage example
"""""""""""""""""""""""

::

   nanopolish eventalign [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``--sam``
     - N
     - NA
     - use to write output in SAM format

   * - ``-w, --window=STR``
     - N
     - NA
     - Compute the consensus for window STR (format : ctg:start_id-end_id)

   * - ``-r, --reads=FILE``
     - Y
     - NA
     - the ONT reads are in fasta FILE

   * - ``-b, --bam=FILE``
     - Y
     - NA
     - the reads aligned to the genome assembly are in bam FILE

   * - ``-g, --genome=FILE``
     - Y
     - NA
     - the genome we are computing a consensus for is in FILE

   * - ``-t, --threads=NUM``
     - N
     - 1
     - use NUM threads

   * - ``--scale-events``
     - N
     - NA
     - scale events to the model, rather than vice-versa

   * - ``--progress``
     - N
     - NA
     - print out a progress message

   * - ``-n``, ``--print-read-names``
     - N
     - NA
     - print read names instead of indexes

   * - ``--summary=FILE``
     - N
     - NA
     - summarize the alignment of each read/strand in FILE

   * - ``--samples``
     - N
     - NA
     - write the raw samples for the event to the tsv output

   * - ``--models-fofn=FILE``
     - N
     - NA
     - read alternative k-mer models from FILE


phase-reads - (experimental)
--------------------

Overview
"""""""""""""""""""""""

Phase reads using heterozygous SNVs with respect to a reference genome 

Input
"""""""""""""""""""""""

    * basecalled reads
    * alignment information
    * assembled genome
    * variants (from nanopolish variants or from other sources eg. Illumina VCF)

Usage example
"""""""""""""""""""""""

::

   nanopolish phase-reads [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa variants.vcf

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * - ``-v``
     - N
     - NA
     - write verbose output

   * - ``-w, --window=STR``
     - N
     - NA
     - Only phase reads in the window STR (format : ctg:start_id-end_id)

   * - ``-r, --reads=FILE``
     - Y
     - NA
     - the ONT reads are in fasta FILE

   * - ``-b, --bam=FILE``
     - Y
     - NA
     - the reads aligned to the genome assembly are in bam FILE

   * - ``-g, --genome=FILE``
     - Y
     - NA
     - the genome we are computing a consensus for is in FILE

   * - ``variants.vcf``
     - Y
     - NA
     - the variants (from nanopolish variants or Illumina in VCF format) to be phased are in FILE

   * - ``-t, --threads=NUM``
     - N
     - 1
     - use NUM threads

    * - ``--progress``
     - N
     - NA
     - print out a progress message

 

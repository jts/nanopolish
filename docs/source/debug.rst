.. _help_us_debug:

Helping us debug nanopolish
===============================

Overview
"""""""""""""""""""""""

Running into errors with nanopolish? To help us debug, we need to be able to reproduce the errors. We can do this by packaging a subset of the files that were used by a nanopolish. We have provided ``scripts/extract_reads_aligned_to_region.py`` and this tutorial to help you do exactly this.

Briefly, this script will:

* extract reads that align to a given region in the draft genome assembly
* rewrite a new BAM, BAI, FASTA file with these reads
* extract the FAST5 files associated with these reads
* save all these files into a tar.gz file

Workflow
"""""""""""""

#. Narrow down a problematic region by running ``nanopolish variants --consensus [...]`` and changing the ``-w`` parameter.
#. Run the ``scripts/extract_reads_aligned_to_region.py``.
#. Send the resulting ``tar.gz`` file to us by hosting either a dropbox or google drive.

.. _creating_example_dataset:

Tutorial - using extraction helper script to create example datsets
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

We extracted a subset of reads for a 2kb region to create the example dataset for the eventalign and consensus tutorial using ``scripts/extract_reads_aligned_to_region.py``. Here is how:

|

Generated basecalled ``--reads`` file:

#. Basecalled reads with albacore: ::

    read_fast5_basecaller.py -c r94_450bps_linear.cfg -t 8 -i /path/to/raw/fast5/files -s /path/to/albacore-2.0.1/output/directory -o fastq 

#. Merged the different albacore fastq outputs: ::

    cat /path/to/albacore-2.0.1/output/directory/workspace/pass/*.fastq > albacore-2.0.1-merged.fastq

#. Converted the merged fastq to fasta format: ::

    paste - - - - < albacore-2.0.1-merged.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > reads.fasta

|

Generated ``--bam`` file with the draft genome assembly (``-g``):

#. Ran canu to create draft genome assembly: ::

    canu \
        -p ecoli -d outdir genomeSize=4.6m \
        -nanopore-raw reads.fasta \ 

#. Index draft assembly: ::

    bwa index ecoli.contigs.fasta
    samtools faidx ecoli.contigs.fasta

#. Aligned reads to draft genome assembly with bwa (v0.7.12): ::

    bwa mem -x ont2d ecoli.contigs.fasta reads.fasta | samtools sort -o reads.sorted.bam -T reads.tmp
    samtools index reads.sorted.bam

|

Selected a ``--window``:

#. Identified the first contig name and chose a random start position: ::

    head -3 ecoli.contigs.fasta

Output: ::

    >tig00000001 len=4376233 reads=23096 covStat=7751.73 gappedBases=no class=contig suggestRepeat=no suggestCircular=no
    AGATGCTTTGAAAGAAACGCAGAATAGATCTCTATGTAATGATATGGAATACTCTGGTATTGTCTGTAAAGATACTAATGGAAAATATTTTGCATCTAAG
    GCAGAAACTGATAATTTAAGAAAGGAGTCATATCCTCTGAAAAGAAAATGTCCCACAGGTACAGATAGAGTTGCTGCTTATCATACTCACGGTGCAGATA
 
As we wanted a 2kb region, we selected a random start position (200000) so our end position was 202000. Therefore our ``--window`` was "tig00000001:200000-202000".

|

Using the files we created, we ran ``scripts/extract_reads_aligned_to_region.py``, please see usage example below.

.. note:: Make sure nanopolish still reproduces the same error on this subset before sending it to us.

Usage example
"""""""""""""""""""""""
::

    python3 extract_reads_aligned_to_region.py \
        --reads reads.fasta \
        --genome ecoli.contigs.fasta \
        --bam reads.sorted.bam \
        --window "tig00000001:200000-202000" \
        --output_prefix ecoli_2kb_region --verbose

.. list-table:: 
   :widths: 25 5 20 50
   :header-rows: 1

   * - Argument name(s)
     - Req.
     - Default value
     - Description

   * - ``-b``, ``--bam``
     - Y
     - NA
     - Sorted bam file created by aligning reads to the draft genome.

   * - ``-g``, ``--genome``
     - Y
     - NA
     - Draft genome assembly

   * - ``-r``, ``--reads``
     - Y
     - NA
     - Fasta, fastq, fasta.gz, or fastq.gz file containing basecalled reads.

   * - ``-w``, ``--window``
     - Y
     - NA
     - Draft genome assembly coordinate string ex. "contig:start-end". It is essential that you wrap the coordinates in quotation marks (\").

   * - ``-o``, ``--output_prefix``
     - N
     - reads_subset
     - Prefix of output tar.gz and log file.

   * - ``-v``, ``--verbose``
     - N
     - False
     - Use for verbose output with info on progress.

Script overview
"""""""""""""""""""""

#. Parse input files
#. Assumes readdb file name from input reads file
#. Validates input
    - checks that input bam, readdb, fasta/q, draft genome assembly, draft genome assembly index file exist, are not empy, and are readable
#. With user input draft genome assembly coordinates, extracts all reads that aligned within these coordinates stores the read_ids. This information can be found from the input BAM.
    - uses pysam.AlignmentFile
    - uses samfile.fetch(region=draft_ga_coords) to get all reads aligned to region
    - if reads map to multiple sections within draft ga it is not added again
#. Parses through the input readdb file to find the FAST5 files associated with each region that aligned to region
    - stores in dictionary region_fast5_files; key = read_id, value = path/to/fast5/file
    - path to fast5 file is currently dependent on the user's directory structure
#. Make a BAM and BAI file for this specific region
    - creates a new BAM file called ``region.bam``
    - with pysam.view we rewrite the new bam with reads that aligned to the region...
    - with pysam.index we create a new BAI file
#. Now to make a new FASTA file with this subset of reads
    - the new fasta file is called ``region.fasta``
    - this first checks what type of sequences file is given { ``fasta``, ``fastq``, ``fasta.gz``, ``fastq.gz`` }
    - then handles based on type of seq file using SeqIO.parse
    - then writes to a new fasta file
#. Let's get to tarring
    - creates a ``tar.gz`` file with the output prefix
    - saves the fast5 files in directory ``output_prefix/fast5_files/``
#. Adds the new fasta, new bam, and new bai file with the subset of reads
#. Adds the draft genome asssembly and associated fai index file
#. Performs a check
    - the number of reads in the new BAM file, new FASTA file, and the number of files in the fast5_files directory should be equal
#. Outputs a ``tar.gz`` and ``log`` file. FIN!

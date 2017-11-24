.. _quickstart:

Quickstart - event align
==================================

The original purpose for Nanopolish was to improve the consensus assembly accuracy for Oxford Nanopore Technology sequencing reads. Here we provide a step-by-step tutorial to help you get started with the nanopolish module: eventalign.

Requirements for tutorial:

* `Nanopolish <installation.html>`_
* `samtool v1.2 <http://samtools.sourceforge.net/>`_
* `bwa mem v0.7.12 <https://github.com/lh3/bwa>`_

Download example dataset
------------------------------------


You can download the example data we will use here: ::

    wget http://s3.climb.ac.uk/nanopolish_tutorial/ecoli_2kb_region.tar.gz
    tar -xvf ecoli_2kb_region.tar.gz
    cd ecoli_2kb_region

Details:

* Sample :    E. coli str. K-12 substr. MG1655
* Instrument : MinION sequencing R9.4 chemistry
* Basecaller : Albacore v2.0.1
* Region: "tig00000001:200000-202000"

This is a subset of reads that aligned to a 2kb region in the E. coli draft assembly.

You should find the following files:

* ``reads.fasta`` : subset of basecalled reads
* ``draft.fa`` : draft genome assembly
* ``draft.fa.fai`` : draft genome assembly index
* ``fast5_files/`` : a directory containing FAST5 files
* ``ecoli_2kb_region.log`` : a log file for how the dataset was created with nanopolish helper script (``scripts/extract_reads_aligned_to_region.py``) 

You will need the reference genome: ::

    wget -O ref.fa ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid225/U00096.ffn


Analysis workflow
------------------


Align the reads with bwa
--------------------------------

First we want to align the reads to the reads: ::

    bwa mem -x ont2d ref.fa reads.fasta | samtools sort -o reads.sorted.bam -T reads.tmp
    samtools index reads.sorted.bam

Align the events with nanopolish eventalign
-----------------------------------------------

Now we align the events to the squiggle data: ::

    nanopolish eventalign \
        --reads reads.fasta \
        --bam reads.sorted.bam \
        --genome ref.fa \
        --sam

This will generate a sam file.

.. _quickstart_eventalign:

Quickstart - how to align events to a reference genome
========================================================

The eventalign module in nanopolish is used to align events or "squiggles" to a reference genome. We (the developers of nanopolish) use this feature extensively when we want to see what the low-level signal information looks like. It helps us model the signal and discover differences in current that might hint at base modifications. Here we provide a step-by-step tutorial to help you get started with the nanopolish eventalign module.

**For more information about eventalign**:

* `Blog post: "Aligning Nanopore Events to a Reference" <http://simpsonlab.github.io/2015/04/08/eventalign/>`_
* `Paper: "A complete bacterial genome assembled de novo using only nanopore sequencing data" <https://www.nature.com/articles/nmeth.3444>`_

**Requirements**:

* `nanopolish <installation.html>`_
* `samtools <http://samtools.sourceforge.net/>`_
* `minimap2 <https://github.com/lh3/minimap2>`_

Download example dataset
------------------------------------

You can download the example dataset we will use here: ::

    wget http://s3.climb.ac.uk/nanopolish_tutorial/ecoli_2kb_region.tar.gz
    tar -xvf ecoli_2kb_region.tar.gz
    cd ecoli_2kb_region

**Details**:

* Sample :    E. coli str. K-12 substr. MG1655
* Instrument : MinION sequencing R9.4 chemistry
* Basecaller : Albacore v2.0.1
* Region: "tig00000001:200000-202000"
* Note: Ligation-mediated PCR amplification performed

This is a subset of reads that aligned to a 2kb region in the E. coli draft assembly. To see how we generated these files please refer to this section: :ref:`creating_example_dataset`. 

You should find the following files:

* ``reads.fasta`` : subset of basecalled reads
* ``fast5_files/`` : a directory containing FAST5 files

You will need the E. coli reference genome: ::

    curl -o ref.fa https://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid225/U00096.ffn

Align the reads with minimap2
--------------------------------

In order to run minimap2 we first need to index the reference genome: ::

    minimap2 -d ref.mmi ref.fa

**Output files**: ``ref.mmi``.

We will need to index the reads as well: ::

    nanopolish index -d fast5_files/ reads.fasta

**Output files**: ``reads.fasta.index``, ``reads.fasta.index.fai``, ``reads.fasta.index.gzi``, and ``reads.fasta.index.readdb``.   

Then we can align the reads to the reference: ::

    minimap2 -ax map-ont -t 8 ref.fa reads.fasta | samtools sort -o reads-ref.sorted.bam -T reads.tmp
    samtools index reads-ref.sorted.bam

**Output files**: ``reads-ref.sorted.bam`` and ``reads-ref.sorted.bam.bai``.

**Checkpoint**: Let's see if the bam file is not truncated. This will check that the beginning of the file contains a valid header, and checks if the EOF is present. This will exit with a non-zero exit code if the conditions were not met: ::

    samtools quickcheck reads-ref.sorted.bam
 
Align the nanopore events to a reference
-----------------------------------------------

Now we are ready to run nanopolish to align the events to the reference genome: ::

    nanopolish eventalign \
        --reads reads.fasta \
        --bam reads-ref.sorted.bam \
        --genome ref.fa \
        --scale-events > reads-ref.eventalign.txt

Assess the eventalign output
-----------------------------------------------

If we take a peek at the first few lines of ``reads-ref.eventalign.txt`` this is what we get: ::

	contig    position  reference_kmer  read_index  strand  event_index  event_level_mean  event_stdv  event_length  model_kmer  model_mean  model_stdv  standardized_level
	gi|545778205|gb|U00096.3|:c514859-514401  3         ATGGAG          0           t       16538        89.82             3.746       0.00100       CTCCAT      92.53       2.49        -0.88
	gi|545778205|gb|U00096.3|:c514859-514401  3         ATGGAG          0           t       16537        88.89             2.185       0.00100       CTCCAT      92.53       2.49        -1.18
	gi|545778205|gb|U00096.3|:c514859-514401  3         ATGGAG          0           t       16536        94.96             2.441       0.00125       CTCCAT      92.53       2.49        0.79
	gi|545778205|gb|U00096.3|:c514859-514401  3         ATGGAG          0           t       16535        81.63             2.760       0.00150       NNNNNN      0.00        0.00        inf
	gi|545778205|gb|U00096.3|:c514859-514401  7         AGTTAA          0           t       16534        78.96             2.278       0.00075       TTAACT      75.55       3.52        0.79
	gi|545778205|gb|U00096.3|:c514859-514401  8         GTTAAT          0           t       16533        98.81             4.001       0.00100       ATTAAC      95.87       3.30        0.72
	gi|545778205|gb|U00096.3|:c514859-514401  9         TTAATG          0           t       16532        96.92             1.506       0.00150       CATTAA      95.43       3.32        0.36
	gi|545778205|gb|U00096.3|:c514859-514401  10        TAATGG          0           t       16531        70.86             0.402       0.00100       CCATTA      68.99       3.70        0.41
	gi|545778205|gb|U00096.3|:c514859-514401  11        AATGGT          0           t       16530        91.24             4.256       0.00175       ACCATT      85.84       2.74        1.60

Example plots
-------------

In `Figure 1 of our methylation detection paper <https://www.nature.com/articles/nmeth.4184>`_ we show a histogram of ``event_level_mean`` for a selection of k-mers to demonstrate how methylation changes the observed current. The data for these figures was generated by eventalign, which we subsequently plotted in R using ggplot2.

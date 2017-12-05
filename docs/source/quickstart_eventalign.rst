.. _quickstart_eventalign:

Quickstart - how to align reads to events
=============================================

The eventalign module in nanopolish is used to align reads to events or "squiggles". As the electrical current signal from MinION has more information, we may be able to more accurately call variants with nanopolish, or compute an improved final consensus sequence. Here we provide a step-by-step tutorial to help you get started with the nanopolish module: eventalign.

**For more information about eventalign**:

* `Blog post: "Aligning Nanopore Events to a Reference" <http://simpsonlab.github.io/2015/04/08/eventalign/>`_
* `Preprint: "A complete bacterial genome assembled de novo using only nanopore sequencing data" <https://www.biorxiv.org/content/early/2015/03/11/015552>`_

**Requirements**:

* `nanopolish v0.8.4 <installation.html>`_
* `samtools v1.2 <http://samtools.sourceforge.net/>`_
* `bwa v0.7.12 <https://github.com/lh3/bwa>`_

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

This is a subset of reads that aligned to a 2kb region in the E. coli draft assembly.

You should find the following files:

* ``reads.fasta`` : subset of basecalled reads
* ``draft.fa`` : draft genome assembly
* ``draft.fa.fai`` : draft genome assembly index
* ``fast5_files/`` : a directory containing FAST5 files
* ``ecoli_2kb_region.log`` : a log file for how the dataset was created with nanopolish helper script (``scripts/extract_reads_aligned_to_region.py``) 

You will need the reference genome: ::

    wget -O ref.fa ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid225/U00096.ffn

Align the reads with bwa
--------------------------------

In order to run bwa with the reference file, we need to index: ::

    bwa index ref.fa

**Output files**: ``ref.fa.sa``, ``ref.fa.amb``, ``ref.fa.ann``, ``ref.fa.pac``, and ``ref.fa.bwt``.

We will need to index the reads as well: ::

    nanopolish index -d fast5_files/ reads.fasta

**Output files**: ``reads.fasta.index``, ``reads.fasta.index.fai``, ``reads.fasta.index.gzi``, and ``reads.fasta.index.readdb``.   

Then we can proceed to aligning the reads to the reference: ::

    bwa mem -x ont2d ref.fa reads.fasta | samtools sort -o reads-ref.sorted.bam -T reads.tmp
    samtools index reads-ref.sorted.bam

**Output files**: ``reads-ref.sorted.bam`` and ``reads-ref.sorted.bam.bai``.

**Checkpoint**: Let's see if the bam file is not truncated. This will check that the beginning of the file contains a valid header, and checks if the EOF is present. This will exit with a non-zero exit code if the conditions were not met: ::

    samtools quickcheck reads-ref.sorted.bam
 
Align the nanopore events to a reference
-----------------------------------------------

Now we are ready to align the events to the squiggle data: ::

    nanopolish eventalign \
        --reads reads.fasta \
        --bam reads-ref.sorted.bam \
        --genome ref.fa > reads-ref.eventalign.txt

Assess the eventalign output
-----------------------------------------------

If we take a peek at the first few lines of ``reads-ref.eventalign.txt`` this is what we get: ::

    contig    position    reference_kmer    read_index    strand    event_index    event_level_mean    event_stdv    event_length    model_kmer    model_mean    model_stdv    standardized_level
    gi|545778205|gb|U00096.3|:c514859-514401    3    ATGGAG    0    t    16538    89.82    3.746    0.00100    CTCCAT    92.53    2.49    -0.88
    gi|545778205|gb|U00096.3|:c514859-514401    3    ATGGAG    0    t    16537    88.89    2.185    0.00100    CTCCAT    92.53    2.49    -1.18
    gi|545778205|gb|U00096.3|:c514859-514401    3    ATGGAG    0    t    16536    94.96    2.441    0.00125    CTCCAT    92.53    2.49    0.79
    gi|545778205|gb|U00096.3|:c514859-514401    3    ATGGAG    0    t    16535    81.63    2.760    0.00150    NNNNNN    0.00    0.00    inf
    gi|545778205|gb|U00096.3|:c514859-514401    7    AGTTAA    0    t    16534    78.96    2.278    0.00075    TTAACT    75.55    3.52    0.79
    gi|545778205|gb|U00096.3|:c514859-514401    8    GTTAAT    0    t    16533    98.81    4.001    0.00100    ATTAAC    95.87    3.30    0.72
    gi|545778205|gb|U00096.3|:c514859-514401    9    TTAATG    0    t    16532    96.92    1.506    0.00150    CATTAA    95.43    3.32    0.36
    gi|545778205|gb|U00096.3|:c514859-514401    10    TAATGG    0    t    16531    70.86    0.402    0.00100    CCATTA    68.99    3.70    0.41
    gi|545778205|gb|U00096.3|:c514859-514401    11    AATGGT    0    t    16530    91.24    4.256    0.00175    ACCATT    85.84    2.74    1.60

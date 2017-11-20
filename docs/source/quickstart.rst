.. _quickstart:

Quickstart
===================

The original purpose for Nanopolish was to improve the consensus assembly accuracy for Oxford Nanopore Technology sequencing reads. Here we provide a step-by-step tutorial to help users get started.

Requirements:

* `Nanopolish <installation.html>`_
* `samtool v1.2 <http://samtools.sourceforge.net/>`_
* `bwa mem v0.7.12 <https://github.com/lh3/bwa>`_

Download example dataset
------------------------------------

You can download an example data set here: ::

    wget http://s3.climb.ac.uk/nanopolish_tutorial/ecoli_2kb_region.tar.gz
	tar -xvf ecoli_2kb_region.tar.gz

Details:

* Sample :	E. coli str. K-12 substr. MG1655
* Instrument : MinION sequencing R9.4 chemistry
* Basecaller : Albacore v2.0.1

You should find the following files:

* ``reads.fasta`` : subset of basecalled reads
* ``draft.fasta`` : draft genome assembly
* ``draft.fasta.fai`` : draft genome assembly index
* ``fast5_files/`` : a directory containing FAST5 files

Overview
-------------------------------
.. image:: _images/nanopolish-workflow.png
   :scale: 50 %
   :alt: nanopolish-tutorial-workflow

   Nanopolish recommended workflow for improving consensus sequence.

Data preprocessing
------------------------------------

Nanopolish needs access to the signal-level data measured by the nanopore sequencer. To begin, we need to create an index ``readdb`` file that links read ids with their signal-level data in the FAST5 files. ::

    nanopolish index -d fast5_files/ reads.fasta

We get the following files: ``reads.fasta.fa.gz``, ``reads.fasta.fa.gz.fai``, ``reads.fasta.fa.gz.gzi``, and ``reads.fasta.fa.gz.readdb``.

Compute the draft genome assembly using CANU
-----------------------------------------------

To create a draft genome assembly we have used CANU. ::

    canu \
		-p ecoli -d outdir genomeSize=4.6m \
		-nanopore-raw albacore-2.0.1-merged.fastq \
		gnuplotTested = true \
		useGrid = false

As this takes a few hours, we have already pre-assembled the data: ``draft.fa``.

Computing a new consensus sequence for a draft assembly
------------------------------------------------------------------------

Now that we have the ``reads.fasta`` indexed with ``nanopolish index``, and have a draft genome assembly ``draft.fa`` we can improve the assembly with nanopolish. The first step is to index the draft genome using BWA INDEX: :: 

    # Index the draft genome
    bwa index draft.fa

Then we align the original non-assembled reads (``reads.fasta``) to the draft assembly (``draft.fa``). ::

    # Align the basecalled reads to the draft sequence
    bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
    samtools index reads.sorted.bam

Now, we use nanopolish to compute the consensus sequence (the genome is polished in 50kb blocks and there will be one output file per block). We'll run this in parallel: ::

    python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fasta -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1

This command will run the consensus algorithm on eight 50kbp segments of the genome at a time, using 4 threads each. Change the ``-P`` and ``--threads`` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual 50kb segments together back into the final assembly: ::

    python nanopolish_merge.py polished.*.fa > polished_genome.fa


Evaluate the assembly
---------------------------------

To analyze how nanopolish performed improving the accuracy we use `MUMmer <https://github.com/mummer4/mummer>`_. MUMmer contains "dnadiff" a script that enables us to see a rreport on alignment statistics. With dnadiff we can compare the two different assemblies. The value that we are interested in is ``AvgIdentity`` which is a measurement of how similar the genome assemblies are to a reference genome.

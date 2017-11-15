.. _quickstart:

Quickstart
===================

Download example dataset
------------------------------------


Data preprocessing
------------------------------------

Nanopolish needs access to the signal-level data measured by the nanopore sequencer. The first step of any nanopolish workflow is to prepare the input data by telling nanopolish where to find the signal files. If you ran Albacore 2.0 on your data you should run nanopolish index on your input reads (``-d`` can be specified more than once if using multiple runs): ::

    # Only run this if you used Albacore 2.0 or later
    nanopolish index -d /path/to/raw_fast5s/ albacore_output.fastq

If you basecalled your reads with Albacore 1.2 or earlier, you should run nanopolish extract on your input reads instead: ::

   # Only run this if you used Albacore 1.2 or later
   nanopolish extract --type template directory/pass/ > reads.fa

Note these two commands are mutually exclusive - you only need to run one of them. You need to decide what command to run depending on the version of Albacore that you used. In the following sections we assume you have preprocessed the data by following the instructions above and that your reads are in a file named ``reads.fa``.

Computing a new consensus sequence for a draft assembly
------------------------------------------------------------------------

The original purpose of nanopolish was to compute an improved consensus sequence for a draft genome assembly produced by a long-read assembly like canu. This section describes how to do this, starting with your draft assembly which should have megabase-sized contigs. ::

    # Index the draft genome
    bwa index draft.fa

    # Align the basecalled reads to the draft sequence
    bwa mem -x ont2d -t 8 draft.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
    samtools index reads.sorted.bam

Now, we use nanopolish to compute the consensus sequence (the genome is polished in 50kb blocks and there will be one output file per block). We'll run this in parallel: ::

    python nanopolish_makerange.py draft.fa | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r reads.fa -b reads.sorted.bam -g draft.fa -t 4 --min-candidate-frequency 0.1

This command will run the consensus algorithm on eight 50kbp segments of the genome at a time, using 4 threads each. Change the ``-P`` and ``--threads`` options as appropriate for the machines you have available.

After all polishing jobs are complete, you can merge the individual 50kb segments together back into the final assembly: ::

    python nanopolish_merge.py polished.*.fa > polished_genome.fa


Calling Methylation
------------------------

nanopolish can use the signal-level information measured by the sequencer to detect 5-mC as described here. Here's how you run it: ::

    # Align the basecalled reads to a reference genome
    bwa mem -x ont2d -t 8 reference.fa reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp -
    samtools index reads.sorted.bam

    # Run the methylation caller
    nanopolish call-methylation -t 8 -r reads.fa -g reference.fa -b reads.sorted.bam > methylation.tsv

The output of call-methylation is a tab-separated file containing per-read log-likelihood ratios (positive values indicate more evidence for 5-mC, negative values indicate more evidence for C). Each line contains the name of the read that covered the CpG site, which allows methylation calls to be phased. We have provided a script to calculate per-site methylation frequencies using call-methylation's output: ::

    python /path/to/nanopolish/scripts/calculate_methylation_frequency -c 2.5 -i methylation.tsv > frequencies.tsv

The output of this script is a tab-seperated file containing the genomic position of the CpG site, the number of reads that covered the site, and the percentage of those reads that were predicted to be methylated. The ``-c 2.5`` option requires the absolute value of the log-likelihood ratio to be at least 2.5 to make a call, otherwise the read will be ignored. This helps reduce calling errors as only sites with sufficient evidence will be included in the calculation.


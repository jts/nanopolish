.. _quickstart_polya:

Quickstart - how to estimate poly(A) tail lengths from nanopore native RNA reads
=================================================================================

Owing to homopolymer effects and the proximity to the sequencing adapters, the polyadenylated tails of reads obtained from nanopore native RNA sequencing are improperly basecalled, making their lengths difficult to quantify. We developed the `polya` subprogram to use an alternative statistical model to estimate these tail lengths.

In this quickstart tutorial, we'll show you how to estimate polyadenylated tail lengths step-by-step, starting from nothing but raw fast5 files. We'll basecall the fast5 files with Oxford Nanopore Technologies' *albacore* basecaller, before aligning the resulting reads with *minimap2*, indexing the files with *nanopolish index*, and finally segmenting the reads and calling the tail lengths with *nanopolish polya*.

We'll be following the steps taken in our `benchmarking analysis <https://github.com/paultsw/polya_analysis>`_ workflow that accompanies `our publication <https://www.biorxiv.org/content/early/2018/11/09/459529>`_. In each step below, we'll generate files in the working directory instead of making subdirectories, except in the case of large datasets.

Software requirements
---------------------
* `nanopolish <https://github.com/jts/nanopolish>`_ `>= v10.2`
* `minimap2 <https://github.com/lh3/minimap2>`_ `>= v2.12`
* `samtools <http://www.htslib.org/>`_ `>= v1.9`
* `albacore >= v2.3.3` (from Oxford Nanopore's private customer portal)

Download raw fast5 data and basecall
------------------------------------
Let's start by downloading a dataset of fast5 files from the European Nucleotide Archive. We'll download a tarball of fast5 files containing reads that are known to have polyadenylated tail lengths of 30nt. ::

    mkdir data && mkdir data/fastqs
    wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA158/ERA1580896/oxfordnanopore_native/30xpolyA.tar.gz -O 30xpolyA.tar.gz && mv 30xpolyA.tar.gz data/
    tar -xzf data/30xpolyA.tar.gz -C data/
    read_fast5_basecaller.py --worker_threads=8 -f FLO-MIN107 -k SQK-RNA001 -q 0 -s data/fastqs -i data/30xpolyA/fast5/pass
    cat data/fastqs/workspace/pass/*.fastq > data/30xpolyA.fastq

In the above, change the value of the `-f` and `-k` arguments based on your flow-cell and sequencing kit, as the basecaller's accuracy is highly dependent upon these settings.

Our directory structure should now look something like this: ::

    (current directory)
    |- data/
    |-- 30xpolyA.fastq
    |-- 30xpolyA.tar.gz
    |-- 30xpolyA/
    |--- fast5/
    |---- pass/
    |----- raw_read1.fast5
    |----- (... more raw fast5's here ...)
    |-- fastqs/
    |--- workspace/
    |---- pass/

Index with nanopolish index
---------------------------
Let's construct an index for our reads with nanopolish's `index` subprogram. This constructs a fast lookup table that tells our program where to find the raw fast5 file for each
basecalled read. ::

    nanopolish index --directory=data/30xpolyA/fast5/pass --sequencing-summary=data/fastqs/sequencing_summary.txt data/30xpolyA.fastq

This should generate a collection of files in the same directory that contains the `30xpolyA.fastq`.

Align with minimap2 and format the BAM file
-------------------------------------------
Now we'll align our basecalled mRNA sequences in the fastq file with our reference. First download a reference fasta: ::

    wget https://raw.githubusercontent.com/paultsw/polya_analysis/master/data/references/enolase_reference.fas
    wget https://raw.githubusercontent.com/paultsw/polya_analysis/master/data/references/enolase_reference.fai
    mv enolase_reference.* data/

Note that our directory structure should now look like this: ::

    (current directory)
    |- data/
    |-- enolase_reference.fas
    |-- enolase_reference.fai
    |-- (... same structure as above ...)

Let's run `minimap2` to align our basecalled reads to this reference and generate a collection of SAM/BAM files with `samtools`: ::

    minimap2 -a -x map-ont data/enolase_reference.fas data/30xpolyA.fastq | samtools view -b - -o data/30xpolyA.bam
    cd data && samtools sort -T tmp -o 30xpolyA.sorted.bam 30xpolyA.bam && samtools index 30xpolyA.sorted.bam && cd ..

Note that we used `-x map-ont`, which is typically for unspliced reads (e.g. genomic DNA) coming from nanopore sequencers. In typical native mRNA sequencing, you should use `-x splice`,
which is a splice-aware setting (and uses a different gap cost in the alignment). We're using `map-ont` due to the fact that our reads come from a control dataset with no splicing.

There should be three more files in the `data` directory now: `30xpolyA.bam`, `30xpolyA.sorted.bam`, and `30xpolyA.sorted.bam.bai`.

Segmentation and tail length estimation with nanopolish polya
-------------------------------------------------------------
Finally, we can run the polyadenylation estimator: ::

    nanopolish polya --threads=8 --reads=data/30xpolyA.fastq --bam=data/30xpolyA.sorted.bam --genome=data/enolase_reference.fas > polya_results.tsv

Set the `--threads` flag to the number of parallel threads you want to use. Generally speaking, a larger number of threads tends to lower the compute time, but there are diminishing
returns to a higher value and performance can actually decrease if your CPU is incapable of supporting your desired number of parallel threads. The best number of threads to use is
highly dependent upon your hardware.

Interpreting the output TSV file
--------------------------------
We'll end this quickstart with a look at the output of the polya program. Let's look at the top five lines of the `polya_results.tsv` file we've just generated: ::

    head -20 polya_results.tsv | column -t
    readname                              contig   position  leader_start  adapter_start  polya_start  transcript_start  read_rate  polya_length  qc_tag
    d6f42b79-90c6-4edd-8c8f-8a7ce0ac6ecb  YHR174W  0         54.0          3446.0         7216.0       8211.0            130.96     38.22         PASS
    453f3f3e-d22f-4d9c-81a6-8576e23390ed  YHR174W  0         228.0         5542.0         10298.0      11046.0           130.96     27.48         PASS
    e02d9858-0c04-4d86-8dba-18a47d9ac005  YHR174W  0         221.0         1812.0         7715.0       8775.0            97.16      29.16         PASS
    b588dee2-2c5b-410c-91e1-fe8140f4f837  YHR174W  0         22.0          8338.0         13432.0      14294.0           130.96     32.43         PASS
    af9dfee2-1711-4083-b109-487b99895e0a  YHR174W  0         889.0         3679.0         6140.0       7168.0            130.96     39.65         PASS
    93f98a86-3b18-48cf-8c4d-15cf277911e2  YHR174W  0         144.0         1464.0         5615.0       6515.0            120.48     30.96         SUFFCLIP
    af9dfee2-1711-4083-b109-487b99895e0a  YHR174W  0         889.0         3679.0         6140.0       7168.0            130.96     39.65         SUFFCLIP
    93f98a86-3b18-48cf-8c4d-15cf277911e2  YHR174W  0         144.0         1464.0         5615.0       6515.0            120.48     30.96         PASS
    ca8d4059-9d82-45ee-aa07-4b8b351618b3  YHR174W  0         1.0           2157.0         4255.0       5862.0            111.56     54.48         PASS
    3493c123-78d4-4f7c-add0-cbb249aef00a  YHR174W  0         78.0          1938.0         4829.0       5491.0            136.91     25.05         PASS
    f5ff1802-3fdd-479a-8888-c72de01bbaea  YHR174W  0         150.0         3476.0         7233.0       7932.0            130.96     25.35         PASS
    bb929728-2ed8-42b0-a5a5-eea4bfd62673  YHR174W  0         91.0          1061.0         6241.0       6910.0            111.56     19.74         PASS
    17cf3fef-1acb-4045-8252-e9c00fedfb7c  YHR174W  0         447.0         6004.0         20058.0      20964.0           100.40     25.17         ADAPTER
    e3e38de6-8a99-4029-a067-261f470517ca  YHR174W  0         41.0          1588.0         4057.0       5303.0            130.96     49.13         PASS
    66f55b56-c22e-4e6d-999e-50687bed6fb7  YHR174W  0         191.0         3160.0         9218.0       10030.0           125.50     28.79         PASS
    56c116d7-9286-4b57-8329-e74928b11b38  YHR174W  0         13.0          1516.0         5845.0       6773.0            130.96     35.30         PASS
    5ca1392c-c48f-4135-85d3-271bd4ee7a13  YHR174W  0         1.0           1976.0         4854.0       5947.0            136.91     44.64         PASS
    66b5a0ef-b8e6-475e-bf20-04b96154a67f  YHR174W  0         98.0          3847.0         7066.0       7925.0            120.48     29.32         PASS
    34bf2187-5816-4744-8e6a-3250b5247e02  YHR174W  0         1.0           2897.0         6885.0       7547.0            125.50     22.54         PASS

Each row corresponds to the output for a given read. The columns have the following interpretation:
* `readname` refers to the unique ID associated to this particular read. This string is also used to look up the corresponding fast5 file, e.g. by looking
  through the `readdb` file generated by `nanopolish index`.
* `contig` refers to the reference sequence that this read aligns to, and is taken from the BAM file.
* `position` is the 5' starting position of the alignment to the reference sequence, and also comes from the BAM file.
* `leader_start`, `adapter_start`, `polya_start`, and `transcript_start` are the indices of the signal segmentation generated by the underlying model within
  `nanopolish`. Briefly, there are four biologically-meaningful regions of the raw sequence of electrical current readings within each fast5 file; these four
  entries denote the starting index of each of these consecutive regions. The indices start from 0 and are oriented in the 3'-to-5' direction (due to the
  inherent orientation of the native RNA nanopore sequencing protocol). A full exposition of this segmentation algorithm is available in the
  `supplementary note<https://www.biorxiv.org/content/biorxiv/suppl/2018/11/09/459529.DC1/459529-2.pdf>`_ to our associated publication.
* `read_rate` is the estimated translocation rate (in units of nucleotides/second) of the read through the pore. The translocation rate is non-uniform during
  the sequencing process of even a single molecule, so this is ultimately a summary statistic of the dynamic, time-varying rate.
* `polya_length` is the estimated polyadenylated tail length, in number of nucleotides. That this value is a float rather than an integer reflects the fact
  that our estimated tail length is the output of an estimator based on the translocation rate.
* `qc_tag` is an additional flag used to indicate the validity of the estimate. Generally speaking, you should only use rows of the output file with this value
  set to `PASS`; all other rows with (e.g.) the `qc_tag` set to `SUFFCLIP`, `ADAPTER`, etc. display signs of irregularity that indicate that we believe the
  estimate to be unreliable. You can easily filter away all rows with the tag set to anything other than `PASS` using a `grep`: ::

    grep 'PASS' polya_results.tsv > polya_results.pass_only.tsv

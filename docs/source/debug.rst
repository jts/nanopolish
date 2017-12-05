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

0. Narrow down a problematic region by running ``nanopolish variants --consensus [...]`` and changing the ``-w`` parameter.
1. Run the ``scripts/extract_reads_aligned_to_region.py``.
2. Send the ``.tar.gz`` to us by hosting either a dropbox or google drive.

Usage example
"""""""""""""""""""""""

::

    python extract_reads_aligned_to_region.py \
        --bam reads.sorted.bam \
        --reads albacore-merged.fastq \
        --genome draft.fa \
		--window "tig00000001:200000-250200" \
        --output_prefix ecoli --verbose

.. list-table:: 
   :widths: 25 5 20 50
   :header-rows: 1

   * - Argument name(s)
     - Req.
     - Default value
     - Description

   * - -b, --bam
     - Y
     - NA
     - Sorted bam file created by aligning reads to the draft genome.

   * - -g, --genome
     - Y
     - NA
     - Draft genome assembly

   * - -r, --reads
     - Y
     - NA
     - Fasta, fastq, fasta.gz, or fastq.gz file containing basecalled reads.

   * - -w, --window
     - Y
     - NA
     - Draft genome assembly coordinate string ex. "contig:start-end". It is essential that you wrap the coordinates in quotation marks (\").

   * - -o, --output_prefix
     - N
     - reads_subset
     - Prefix of output tar.gz and log file.

   * - -v, --verbose
     - N
     - False
     - Use for verbose output with info on progress.

Script overview
"""""""""""""""""""""

0. Parse input files
1. Assumes readdb file name from input reads file
2. Validates input
    - checks that input bam, readdb, fasta/q, draft genome assembly, draft genome assembly index file exist, are not empy, and are readable
3. With user input draft genome assembly coordinates, extracts all reads that aligned within these coordinates
   stores the read_ids. This information can be found from the input BAM.
    - uses pysam.AlignmentFile
    - uses samfile.fetch(region=draft_ga_coords) to get all reads aligned to region
    - if reads map to multiple sections within draft ga it is not added again
4. Parses through the input readdb file to find the FAST5 files associated with each region that aligned to region
   - stores in dictionary region_fast5_files; key = read_id, value = path/to/fast5/file
   - path to fast5 file is currently dependent on the user's directory structure
5. Make a BAM and BAI file for this specific region
   - creates a new BAM file called region.bam
   - with pysam.view we rewrite the new bam with reads that aligned to the region...
   - with pysam.index we create a new BAI file
6. Now to make a new FASTA file with this subset of reads
   - the new fasta file is called region.fasta
   - this first checks what type of sequences file is given { fasta, fastq, fasta.gz, fastq.gz }
   - then handles based on type of seq file using SeqIO.parse
   - then writes to a new fasta file
7. Let's get to tarring
   - creates a tar.gz file with the output prefix
   - saves the fast5 files in directory output_prefix/fast5_files/
8. Adds the new fasta, new bam, and new bai file with the subset of reads
9. Adds the draft genome asssembly and associated fai index file
10. Performs a check
   - the number of reads in the new BAM file, new FASTA file, and the number of files in the fast5_files directory should be equal
11. Outputs a log file
FIN!

import pysam
import sys
import re
import subprocess
import os
from Clustal import *
from collections import defaultdict
from Bio import AlignIO

# reverse complement a sequence
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def cleanup_temp(fn):
    print 'not deleting', fn
    #os.remove(fn)

# parse an LAshow read index string into a numeric ID
# the IDs are 1-based and contain commas
def lashow_idstr2idx(s):
    return int(s.replace(',', '')) - 1

# remove non-numeric characters from a string
def remove_nonnumeric(s):
    return re.sub("[^0-9]", "", s)

# parse an LAshow output file and build a map from a read index
# to the sequences that align to it
def parse_lashow(fn):
    fh = open(fn, 'r')
    
    out = defaultdict(list)

    for line in fh:

        fields = line.split()
        if len(fields) != 18:
            continue
       
        id1 = lashow_idstr2idx(fields[0])
        id2 = lashow_idstr2idx(fields[1])
        strand = fields[2]

        # 
        s = int(remove_nonnumeric(fields[8]))
        e = int(remove_nonnumeric(fields[9]))
    
        out[id1].append((id2, strand, s, e))
    return out

# return true if the pair of intervals intersect
def do_regions_overlap(s1, e1, s2, e2):
    i1 = s1 >= s2 and s1 <= e2
    i2 = s2 >= s1 and s2 <= e1
    return i1 or i2

# write a fasta file for input into POA
def write_poa_input(input_list, segment_idx):
    fn = "poa.input.%d.fa" % (segment_idx)
    fh = open(fn, "w")

    n_reads = 0
    for poa_input in input_list:
        poa_read = poa_input.read
        if not poa_read.is_base:
            out_id = "%d:%c:%d-%d" % (poa_read.read_id, poa_read.strand, poa_read.start, poa_read.stop)
        else:
            out_id = "%d:%s" % (segment_idx, "poabaseread")
        fh.write(">%s\n%s\n" % (out_id, poa_input.sequence))
        n_reads += 1

    fh.close()
    return (fn, n_reads)

#
def run_poa_and_consensus(input_list, segment_idx):
    (in_fn, n_reads) = write_poa_input(input_list, segment_idx)
    out_fn = "clustal-%d.out" % (segment_idx)
    cmd = "poa -read_fasta %s -clustal %s -hb poa-blosum80.mat" % (in_fn, out_fn)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    
    #consensus =  clustal2consensus(out_fn)

    cleanup_temp(in_fn)
    cleanup_temp(out_fn)
    return ("", n_reads)

# Args
assembly_fn = sys.argv[1]
bam_fn = sys.argv[2]

# Build an index from fast5 file names to indices in the map file
fast5_fofn_fn = 'r73.map.fofn'
fast5_fofn_fh = open(fast5_fofn_fn)

f5_index_map = dict()
for (idx, file_path) in enumerate(fast5_fofn_fh):
    # Strip the leading directories from the name
    f5_fn = os.path.basename(file_path.rstrip())
    f5_index_map[f5_fn] = idx

# Open assembly file
assembly = pysam.Fastafile(assembly_fn)

# Open bam file
bamfile = pysam.AlignmentFile(bam_fn, "rb")

# Iterate over contigs
SEGMENT_LENGTH = 10000
print '!!!!!!!!!! SKIPPING START !!!!!!!!!!'

segment_id = 0
for (contig_idx, contig_name) in enumerate(bamfile.references):
    contig_length = bamfile.lengths[contig_idx]
    print 'Computing consensus for', contig_name, contig_length, 'bp'

    for segment_start in xrange(20000, contig_length, SEGMENT_LENGTH):
        segment_end = segment_start + SEGMENT_LENGTH
        if segment_end > contig_length:
            segment_end = contig_length
    
        segment_sequence = assembly.fetch(contig_name, segment_start, segment_end)
        input_list = list()

        # Make a sequence for the assembly segment
        poa_read = PoaRead(-1, 'n', 0, -1, True)
        poa_input = PoaInput(poa_read, segment_sequence)
        input_list.append(poa_input)

        # Iterate over the bases in the segment and extract the reads aligned here
        for r in bamfile.fetch(contig_name, segment_start, segment_end):
            
            print r.query_name, r.query_alignment_start, r.query_alignment_end, r.reference_start, r.reference_end

            # Clip the reads coordinates to this segment
            first_read_base_in_segment = -1
            last_read_base_in_segment = -1
            for (read_idx, reference_idx) in enumerate(r.get_reference_positions(full_length=True)):

                if reference_idx >= segment_start and first_read_base_in_segment == -1:
                    first_read_base_in_segment = read_idx
                
                if reference_idx is not None and reference_idx < segment_end:
                    last_read_base_in_segment = read_idx
            
            #
            print "\t", r.query_name, first_read_base_in_segment, last_read_base_in_segment
        
            # read sub sequence
            sequence = r.query_sequence[first_read_base_in_segment:last_read_base_in_segment]

            orientation = 'n'
            if r.is_reverse:
                orientation = 'c'
 
            # Get the index of the read in the file map
            # the consensus caller needs this and POA truncates long names
            f5_fn = os.path.basename(r.query_name)
            if f5_fn not in f5_index_map:
                print 'Cannot find filename for', f5_fn

            read_idx = f5_index_map[f5_fn]
            poa_read = PoaRead(read_idx, orientation, first_read_base_in_segment, last_read_base_in_segment, False)
            poa_input = PoaInput(poa_read, sequence)
            input_list.append(poa_input)

        run_poa_and_consensus(input_list, segment_id)
        segment_id += 1
        sys.exit(1)

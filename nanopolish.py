import pysam
import sys
import re
import subprocess
import os
import random
import argparse
from Clustal import *
from HmmCons import *
from collections import defaultdict
from Bio import AlignIO
from Bio import pairwise2

# reverse complement a sequence
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def cleanup_temp(fn):
    #print "Cleanup turned off"
    os.remove(fn)

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
def run_poa_and_consensus(input_list, segment_idx, num_threads):
    (in_fn, n_reads) = write_poa_input(input_list, segment_idx)
    out_fn = "clustal-%d.out" % (segment_idx)
    cmd = "poa -read_fasta %s -clustal %s -hb poa-blosum80.mat > /dev/null" % (in_fn, out_fn)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    
    consensus = call_consensus_for_file(out_fn, num_threads)
    
    # fake a consensus sequence for debugging
    #consensus = input_list[0].sequence

    cleanup_temp(in_fn)
    cleanup_temp(out_fn)
    return (consensus, n_reads)

def merge_into_consensus(consensus, incoming, overlap_length):

    # if first segment, no overlapping needs to be done
    if consensus == "":
        return incoming

    or_con = consensus[-overlap_length:]
    or_inc = incoming[0:overlap_length]

    # These parameters are designed to give us the region of highest similarity
    # between the two sequences
    alignments = pairwise2.align.localms(or_con, or_inc, 2, -10, -10, -10)

    best = alignments[0]
    aln_con, aln_inc, score, begin, end = best

    # Calculate the point to merge at
    m_con = 0
    m_inc = 0

    for i in xrange(0, begin):
        a = aln_con[i]
        b = aln_inc[i]

        if a != '-':
            m_con += 1
        if b != '-':
            m_inc += 1

    m_con += len(consensus) - overlap_length
    merged = consensus[0:m_con] + incoming[m_inc:]
    
    #merged = or_con[0:m_con] + or_inc[m_inc:]
    debug_segment_merge = False
    if debug_segment_merge:
        print 'OR', or_con
        print 'OR', or_inc

        print 'Before trim'
        print aln_con
        print aln_inc

        print 'After trim'
        print aln_con[begin:end]
        print aln_inc[begin:end]
        print score, begin, end, m_con, m_inc
        print merged

    return merged

# Args
parser = argparse.ArgumentParser()
parser.add_argument("--assembly", help="the filename of the assembly to polish")
parser.add_argument("--bam", help="the name of the bam containing reads mapped to the assembly")
parser.add_argument("--threads", help="the number of compute threads to use")
parser.add_argument("--contig", help="only compute the consensus for the given contig")
parser.add_argument("--segments", help="only compute the consensus for the given range of segments")
args = parser.parse_args()

# Open assembly file
assembly = pysam.Fastafile(args.assembly)

# Open bam file
bamfile = pysam.AlignmentFile(args.bam, "rb")

# Open partial output file
#final_consensus_fh = open("nanopolish_final.fa", "w")
contig_out_name = "all"
if args.contig is not None:
    contig_out_name = args.contig

segment_out_name = "all"
if args.segments is not None:
    segment_out_name = args.segments

segment_consensus_fh = open("nanopolish_contig_%s_segment_%s.fa" % (contig_out_name, segment_out_name), "w")

# Iterate over contigs
SEGMENT_LENGTH = 10000
SEGMENT_OVERLAP = 200

for (contig_idx, contig_name) in enumerate(bamfile.references):

    if args.contig != None and args.contig != contig_name:
        continue

    contig_length = bamfile.lengths[contig_idx]
    print 'Computing consensus for', contig_name, contig_length, 'bp'
    contig_consensus = ""

    num_segments = int(contig_length / SEGMENT_LENGTH) + 1
    segments_to_process = list()
    if args.segments == None:
        # process entire contig
        segments_to_process = range(0, num_segments)
    else:
        # process range of segments
        (start, end) = [ int(x) for x in args.segments.split(':') ]
        if end > num_segments:
            end = num_segments

        segments_to_process = range(start, end)

    for segment_id in segments_to_process:
        segment_start = segment_id * SEGMENT_LENGTH
        segment_end = segment_start + SEGMENT_LENGTH + SEGMENT_OVERLAP

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
            
            # Clip the reads coordinates to this segment
            first_read_base_in_segment = -1
            last_read_base_in_segment = -1
            for (read_idx, reference_idx) in enumerate(r.get_reference_positions(full_length=True)):

                if reference_idx >= segment_start and first_read_base_in_segment == -1:
                    first_read_base_in_segment = read_idx
                
                if reference_idx is not None and reference_idx < segment_end:
                    last_read_base_in_segment = read_idx
            
            # read sub sequence
            sequence = r.query_sequence[first_read_base_in_segment:last_read_base_in_segment]

            orientation = 'n'
            if r.is_reverse:
                orientation = 'c'
            
            # Parse the read index out of the query name
            read_idx = int(r.query_name.split('.')[0])
            poa_read = PoaRead(read_idx, orientation, first_read_base_in_segment, last_read_base_in_segment, False)
            poa_input = PoaInput(poa_read, sequence)
            input_list.append(poa_input)

        sys.stderr.write("\tcalculating consensus for segment %d:%d\n" % (segment_start, segment_end))

        (segment_consensus, num_reads) = run_poa_and_consensus(input_list, segment_id, args.threads)

        # Write segment consensus to file
        segment_consensus_fh.write(">%s:%s\n%s\n" % (contig_name, segment_id, segment_consensus))
        #contig_consensus = merge_into_consensus(contig_consensus, segment_consensus, SEGMENT_OVERLAP)
        segment_id += 1
    
    # Write final consensus to file
    #final_consensus_fh.write(">%s\n%s\n" % (contig_name, contig_consensus))

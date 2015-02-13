import sys
import glob
from Bio import SeqIO
from Bio import pairwise2

def merge_into_consensus(consensus, incoming, overlap_length):

    # if first segment, no overlapping needs to be done
    if consensus == "":
        return incoming

    or_con = consensus[-overlap_length:]
    or_inc = incoming[0:overlap_length]

    # These parameters are designed to give us the region of highest similarity
    # between the two sequences
    alignments = pairwise2.align.globalms(or_con, or_inc, 2, -10, -10, -3)

    best = alignments[0]
    aln_con, aln_inc, score, begin, end = best

    # We merge at the midpoint between the two aligned segments
    m_con = 0
    m_inc = 0

    assert(len(aln_con) == len(aln_inc))

    for i in xrange(0, len(aln_con) / 2):
        a = aln_con[i]
        b = aln_inc[i]

        if a != '-':
            m_con += 1
        if b != '-':
            m_inc += 1

    
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

        print 'Merging:'
        print or_con[0:m_con]
        print or_inc[m_inc:]

    m_con += len(consensus) - overlap_length
    merged = consensus[0:m_con] + incoming[m_inc:]

    return merged

# Make placeholder segments using the original assembly as a guide
original_assembly = sys.argv[1]
recs = [ (rec.name, len(rec.seq)) for rec in SeqIO.parse(open(original_assembly), "fasta")]

# Do not change, must match nanopolish segment lengths
SEGMENT_LENGTH = 10000
OVERLAP_LENGTH = 200

segments_by_name = dict()
for name, length in recs:

    n_segments = (length / SEGMENT_LENGTH) + 1
    segments_by_name[name] = [""] * n_segments

for fn in sys.argv[2:]:
    for rec in SeqIO.parse(open(fn), "fasta"):
        (contig, segment) = rec.name.split(":")
        segments_by_name[contig][int(segment)] = str(rec.seq)

# Assemble while making sure every segment is present
for contig_name in sorted(segments_by_name.keys()):
    assembly = ""
    for (segment_id, sequence) in enumerate(segments_by_name[contig_name]):
        
        if sequence is "":
            print "ERROR, segment %d of contig %s is missing" % (segment_id, contig_name)
            sys.exit(1)

        sys.stderr.write('Merging %s %d\n' % (contig_name, segment_id))

        assembly = merge_into_consensus(assembly, sequence, OVERLAP_LENGTH)
    
    # Write final assembly
    print(">%s\n%s\n" % (contig_name, assembly))

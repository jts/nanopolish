import sys
from Bio import SeqIO

recs = [ (rec.name, len(rec.seq)) for rec in SeqIO.parse(open(sys.argv[1]), "fasta")]

# Do not change, must match nanopolish segment lengths
SEGMENT_LENGTH = 10000

# Ok to change this
SEGMENTS_PER_BATCH = 10

for name, length in recs:
    n_segments = (length / SEGMENT_LENGTH) + 1

    for n in xrange(0, n_segments, SEGMENTS_PER_BATCH):
        if ( n + SEGMENTS_PER_BATCH) > n_segments:
            print "%s:%d-%d" % (name, n, n_segments)
        else:
            print "%s:%d-%d" % (name, n, n + SEGMENTS_PER_BATCH)

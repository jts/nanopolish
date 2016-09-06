import sys
from Bio import SeqIO

recs = [ (rec.name, len(rec.seq)) for rec in SeqIO.parse(open(sys.argv[1]), "fasta")]

SEGMENT_LENGTH = 50000
OVERLAP_LENGTH = 200

for name, length in recs:
    n_segments = (length / SEGMENT_LENGTH) + 1

    for n in xrange(0, length, SEGMENT_LENGTH):
        if ( n + SEGMENT_LENGTH) > length:
            print "%s:%d-%d" % (name, n, length - 1)
        else:
            print "%s:%d-%d" % (name, n, n + SEGMENT_LENGTH + OVERLAP_LENGTH)

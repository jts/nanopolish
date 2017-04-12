import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Partition a genome into a set of overlapping segments')
parser.add_argument('--segment-length', type=int, default=50000)
parser.add_argument('--overlap-length', type=int, default=200)
args, extra = parser.parse_known_args()
if len(extra) != 1:
    sys.stderr.write("Error: a genome file is expected\n")
filename = extra[0]

recs = [ (rec.name, len(rec.seq)) for rec in SeqIO.parse(open(filename), "fasta")]

SEGMENT_LENGTH = args.segment_length
OVERLAP_LENGTH = args.overlap_length

for name, length in recs:
    n_segments = (length / SEGMENT_LENGTH) + 1

    for n in xrange(0, length, SEGMENT_LENGTH):
        if ( n + SEGMENT_LENGTH) > length:
            print("%s:%d-%d" % (name, n, length - 1))
        else:
            print("%s:%d-%d" % (name, n, n + SEGMENT_LENGTH + OVERLAP_LENGTH))

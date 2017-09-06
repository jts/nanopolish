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
MIN_SEGMENT_LENGTH = 5 * OVERLAP_LENGTH

for name, length in recs:
    n_segments = (length / SEGMENT_LENGTH) + 1

    start = 0
    while start < length:
        end = start + SEGMENT_LENGTH

        # If this segment will end near the end of the contig, extend it to end
        if length - end < MIN_SEGMENT_LENGTH:
            print("%s:%d-%d" % (name, start, length - 1))
            start = length
        else:
            print("%s:%d-%d" % (name, start, end + OVERLAP_LENGTH))
            start = end

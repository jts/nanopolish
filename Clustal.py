import sys
from collections import namedtuple
from NanoUtils import *
from SquiggleRead import *
from Bio import AlignIO

PoaRead = namedtuple('PoaRead', ['index', 'strand', 'start', 'stop'])

# Parse a read id to demangle the meta data encoded within it
def unpack_poa_id(rid):
    a = rid.split(':')
    if a[1] == 'poabaseread':
        
        return PoaRead(int(a[0]), 'n', 0, -1)
    else:
        c = a[2].split('-')
    return PoaRead(int(a[0]), a[1], int(c[0]), int(c[1]))

class Clustal:

    def __init__(self, filename):
        self.alignment = AlignIO.read(filename, "clustal")

        # Get the read names out of the alignment rows
        read_ids = self.get_read_ids()

        # Unpack the metadata from the reads
        self.reads = []
        for rid in read_ids:
            self.reads.append(unpack_poa_id(rid))

    # Return the index of the row giving the first consensus sequence
    def get_consensus_row(self):
        for (i, record) in enumerate(self.alignment):
            if record.id == 'CONSENS0':
                return i
        return -1

    # Return the indices of non-consensus rows
    def get_read_rows(self):
        out = []
        for (i, record) in enumerate(self.alignment):
            if record.id.find('CONSEN') == -1:
                out.append(i)
        return out

    # Returns the IDs of the records that are reads
    def get_read_ids(self):
        indices = self.get_read_rows()
        return [ self.alignment[x].id for x in indices ]

    # subset the multiple alignment to only show the i-th and j-th rows
    def show_pairwise(self, i, j):
        ri = self.alignment[i]
        rj = self.alignment[j]

        match = 0
        mismatch = 0
        gap = 0

        # iterate over columns
        oi = []
        oj = []
        for (a, b) in zip(ri, rj):

            # do nothing when both are a gap
            if a == '-' and b == '-':
                continue

            if a == b:
                match += 1
            elif a == '-' or b == '-':
                gap += 1
            else:
                mismatch += 1

            oi.append(a)
            oj.append(b)

        print ''.join(oi)
        print ''.join(oj)
        print match,mismatch,gap,(float(match) / (match + mismatch + gap))

    # Print the indices of 5mers that are shown at the same aligned
    # position of all rows
    def print_consensus_5mers(self):
        n_rows = len(self.alignment)
        n_cols = len(self.alignment[0]) # should be same for all rows
        n_mers = n_cols - 4

        # An array containing the number of bases seen for each row so far
        base_indices = [0] * n_rows

        for i in xrange(0, n_mers):
            k0 = str(self.alignment[0][i:i+5].seq)

            all_match = True
            for j in xrange(1, n_rows):
                kj = str(self.alignment[j][i:i+5].seq)
                if k0 != kj:
                    all_match = False
                    break
            

            if all_match:
                # k-mer shared between all rows should never contain a gap
                if k0.find('-') != -1:
                    print 'Error, gapped k-mer match:', k0, i
                    sys.exit(1)

                out = []
                out.append(k0)
                for j in xrange(0, n_rows):
                    out.append(str(base_indices[j]))
                print '\t'.join(out)

            # Update base indices
            for j in xrange(0, n_rows):
                if self.alignment[j][i] != '-':
                    base_indices[j] += 1

    # Get k-mers from the consensus sequence that are supported
    # by at least min_count reads
    def get_supported_consensus_kmers(self, min_count):

        row_idx = self.get_consensus_row()
        read_row_indices = self.get_read_rows()
        print read_row_indices
        n_rows = len(self.alignment)
        n_cols = len(self.alignment[0]) # should be same for all rows
        n_mers = n_cols - 4

        for i in xrange(0, n_mers):

            # Skip columns that are consensus gaps
            if self.alignment[row_idx][i] == '-':
                continue

            mer = self.get_kmer(row_idx, i, 5)

            # Count the number of reads that have the same pattern in these columns
            support = 0
            for j in read_row_indices:
               if str(self.alignment[j][i:stop].seq) == mer:
                   support += 1
            print i, mer, support

    # Get at least k bases starting from a given row and column 
    def get_kmer(self, row, col, k):
        # Extract bases from this column until
        # we have 5 bases. The string can be longer
        # than 5 due to the precense of gaps
        bc = 0
        stop = col
        while bc < 5 and stop < len(self.alignment[row]):
            if self.alignment[row][stop] != '-':
                bc += 1
            stop += 1
        return str(self.alignment[row][col:stop].seq)

    # Count the number of bases seen up to and including the given column
    def count_bases_to_column(self, row, col):
        bc = 0
        for i in xrange(0, col):
            bc += (self.alignment[row][i] != '-')
        return bc

    # Check whether each read matches the given k-mer at the specified column
    # Returns a list of positions in the original reads of the matched k-mer
    def match_column_kmer(self, col, kmer):
        out = []
        for ri in self.get_read_rows():

            match = True
            for i in xrange(0, len(kmer)):
                if self.alignment[ri][col + i] != kmer[i]:
                    match = False
                    break
            if match:
                # SLOW
                base_index = self.count_bases_to_column(ri, col) + self.reads[ri].start
                out.append((ri, base_index))
            else:
                out.append((ri, -1))
        return out

    # Compare a pair of consensus rows
    def compare_consensus_rows(self):
        con0_idx = self.get_consensus_row()
        con1_idx = con0_idx + 1

        n_rows = len(self.alignment)
        n_cols = len(self.alignment[0]) # should be same for all rows
        n_mers = n_cols - 4
        
        matches = []
        branches = []

        for i in xrange(0, n_mers):
            k0 = self.get_kmer(con0_idx, i, 5)
            k1 = self.get_kmer(con1_idx, i, 5)

            if k0 == k1:
                matches.append(k0)
            else:
                branches.append((i - 1, matches[-1], k0, k1))

        return branches

if __name__ == '__main__':

    ca_fn = 'clustal-32.out'
    fast5_fofn_fn = 'r73.map.fofn'
    fast5_fofn_fh = open(fast5_fofn_fn)

    files = []
    for l in fast5_fofn_fh:
        files.append(l.rstrip())
     
    ca = Clustal(ca_fn)
    read_id_to_squiggle = {}

    squiggle_reads = []
    for read in ca.reads:
        print read, files[read.index]
        squiggle_reads.append(SquiggleRead(files[read.index]))

    #ca.show_pairwise(0, 1)
    #ca.show_pairwise(0, 2)
    #ca.print_consensus_5mers()
    #ca.get_supported_consensus_kmers(3)
    branches = ca.compare_consensus_rows()
    b = branches[0]
    print b

    matches = ca.match_column_kmer(b[0], b[1])
    for (idx, k_idx) in matches:

        if k_idx == -1:
            continue

        read = ca.reads[idx]
        
        sr = squiggle_reads[idx]
        
        # Setup branch k-mers
        x = branches[0][1]
        y = branches[0][2]
        z = branches[0][3]

        # If these read was used in reverse-complement orientation
        # to build the POA, we need to reverse the coordinates and k-mers below
        if read.strand == 'c':
            x = revcomp(x)
            y = revcomp(y)
            z = revcomp(z)
            k_idx = sr.flip_k_idx_strand(k_idx, 5)
            (branch_te, branch_ce) = sr.get_event_before_2D_kmer(k_idx)
        else:
            (branch_te, branch_ce) = sr.get_event_after_2D_kmer(k_idx)

        print 'Computing branches:', read.index, read.strand, x, sr.get_2D_kmer_at(k_idx, 5)

        print 'Matched:'
        matched_events = sr.get_events_for_2D_kmer(k_idx)
        for (te, ce) in matched_events:
            print sr.format_kmer_event_pair(te, ce, x)

        print 'Next:'
        print sr.format_kmer_event_pair(branch_te, branch_ce, y)
        print sr.format_kmer_event_pair(branch_te, branch_ce, z)

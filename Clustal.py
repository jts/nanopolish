import sys
from collections import namedtuple, defaultdict
from NanoUtils import *
from SquiggleRead import *
from Bio import AlignIO

PoaRead = namedtuple('PoaRead', ['read_id', 'strand', 'start', 'stop', 'is_base'])
PoaInput = namedtuple('PoaInput', ['read', 'sequence'])

# Parse a read id to demangle the meta data encoded within it
def unpack_poa_id(rid):
    a = rid.split(':')
    if a[1] == 'poabaseread':
        return PoaRead(int(a[0]), 'n', 0, -1, 1)
    else:
        c = a[2].split('-')
    return PoaRead(int(a[0]), a[1], int(c[0]), int(c[1]), 0)

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

    # Get at least k bases starting from a given row and column 
    def get_kmer(self, row, col, k):
        # Extract bases from this column until
        # we have 5 bases. The string can be longer
        # than 5 due to the precense of gaps
        bc = 0
        stop = col
        while bc < k and stop < len(self.alignment[row]):
            if self.alignment[row][stop] != '-':
                bc += 1
            stop += 1
        return str(self.alignment[row][col:stop].seq).replace('-', '')

    # Get the sequence in the given row spanning the columns
    # up to and including the kmer starting in the last column
    def get_sequence_plus_k(self, row, start_col, end_col, K):
        
        # Get the sequence up to the last column
        # and append the following k-1 bases
        seq = str(self.alignment[row][start_col:end_col].seq).replace('-', '') + self.get_kmer(row, end_col, K)
        return seq

    # Convert all sequences to upper case strings
    # and remove Ns
    def preprocess(self):
        for row in self.alignment:
            row.seq = str(row.seq).upper()
            ns = list()
            for i in xrange(0, len(row.seq)):
                if row.seq[i] == 'N' or row.seq[i] == 'n':
                    ns.append(random_base())
                else:
                    ns.append(row.seq[i])
            row.seq = ''.join(ns)
    
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

    # Returns the indices of reads that have a base (not a gap)
    # in the specified column
    def get_reads_with_bases_in_columns(self, i, j):
        out = []
        for ri in self.get_read_rows():
            if self.alignment[ri][i] != '-' and self.alignment[ri][j] != '-':
                out.append(ri)
        return out

    # Generate the possible consensus sequences from
    # column i to column j using the reads in read_indices
    def generate_possible_consensus(self, read_indices, i, j):
        out = [""]
        while i <= j:
            symbols = defaultdict(int)
            for ri in read_indices:
                symbols[self.alignment[ri][i]] += 1

            new_out = []

            # Append every symbol to every previous string
            for (base, count) in symbols.iteritems():
                #if count <= 1:
                #    continue

                for seq in out:
                    if base != '-':
                        new_out.append(seq + base)
                    else:
                        new_out.append(seq)
            out = new_out
            i += 1
        return out

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
        #squiggle_reads.append(SquiggleRead(files[read.index]))

    #ca.show_pairwise(0, 1)
    #ca.show_pairwise(0, 2)
    #ca.print_consensus_5mers()

    col_start = 0
    col_stop = 10
    read_indices = ca.get_reads_with_bases_in_columns(col_start, col_stop)
    consensus_set = ca.generate_possible_consensus(read_indices, col_start, col_stop)
    print consensus_set

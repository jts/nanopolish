from SquiggleRead import *
from NanoUtils import *
from Clustal import *
from ctypes import *
from collections import defaultdict, namedtuple

StrandAnchor = namedtuple('StrandAnchor', ['idstr', 'row_id', 'sr_idx', 'strand', 'event_idx', 'rc', 'diff'])

class Anchor:
    def __init__(self, column, cons_kidx):
        self.column = column
        self.cons_kidx = cons_kidx
        self.strands = list()
        self.n_outliers = 0
        self.n_samples = 0
        self.sum_diff = 0

    def add_strand(self, sa):
        self.strands.append(sa)

# Convenience type for working with ctypes
c_double_p = POINTER(c_double)

class CPoreModelInterface(Structure):
    _fields_ = [('n_states', c_int),
                ('pore_model_level_mean', c_double_p),
                ('pore_model_level_stdv', c_double_p),
                ('pore_model_sd_mean', c_double_p),
                ('pore_model_sd_stdv', c_double_p),
                ('pore_model_scale', c_double),
                ('pore_model_shift', c_double),
                ('pore_model_drift', c_double),
                ('pore_model_var', c_double)]

class CEventSequenceInterface(Structure):
    _fields_ = [('n_events', c_int),
                ('level', c_double_p),
                ('stdv', c_double_p),
                ('time', c_double_p)]

class CSquiggleReadInterface(Structure):
    _fields_ = [('pore_model', CPoreModelInterface * 2),
                ('events', CEventSequenceInterface * 2)]

class CReadStateInterface(Structure):
    _fields_ = [('read_idx', c_int),
                ('event_start_idx', c_int),
                ('event_stop_idx', c_int),
                ('strand', c_int),
                ('stride', c_int),
                ('rc', c_int)]

class CReadAnchorInterface(Structure):
    _fields_ = [('event_idx', c_int),
                ('rc', c_int)]

# Generate all the paths of length (k+1)
# starting from the given k-mer
def generate_kmer_paths(kmer):
    k = len(kmer)
    all_kmers = generate_mers(dna_alphabet, k)

    out = []
    for a in all_kmers:
        path_str = kmer + a
        out.append(str2kmers(kmer + a, k))
    return out

# get the next observed event in the range
def get_next_event(sr, start, stop, stride, strand):
    for i in xrange(start, stop, stride):

        events = sr.event_map['2D'][i]
        if len(events) == 0:
            continue

        ei = -1

        if strand == 't':
            ei = events[0][0]
        else:
            ei = events[0][1]
        
        if ei != -1:
            return ei
    return -1

# Get the closest event to the indicated
# 2D kmer for the given squiggle read
def get_closest_event_to(sr, kidx, strand, kmer):
    e_level = sr.get_expected_level(kmer, strand)
    
    first_event_before = -1
    first_event_after = -1
    m = 1000
    
    # before
    stop = kidx - m
    if stop < 0:
        stop = -1 # end is exclusive
    
    first_event_before = get_next_event(sr, kidx, stop, -1, strand)

    # after
    stop = kidx + m
    if stop >= len(sr.event_map['2D']):
        stop = len(sr.event_map['2D'])

    first_event_after = get_next_event(sr, kidx, stop, 1, strand)

    if first_event_before is None or first_event_after is None:
        return -1

    # evaluate levels
    b_level = sr.get_drift_corrected_event_level(first_event_before, strand)
    a_level = sr.get_drift_corrected_event_level(first_event_after, strand)
    #print 'SEARCH', first_event_before, first_event_after, b_level, a_level, e_level

    if abs(b_level - e_level) < abs(a_level - e_level):
        return first_event_before
    else:
        return first_event_after


def calculate_diff_to_expected(sr, kidx, strand, kmer):
    
    ei = get_closest_event_to(sr, kidx, strand, kmer)
    if ei == -1:
        return (float("-inf"), ei)

    # Observed event   
    level = sr.get_drift_corrected_event_level(ei, strand)

    k_level = sr.get_expected_level(kmer, strand)
    k_sd = sr.get_expected_sd(kmer, strand)
    res = (level - k_level) / k_sd
    return (res, ei)

# For the pair of input Anchors, make a list
# of StrandAnchors thare are anchored at both inputs
def pair_anchor_strands(a1, a2):

    out = list()

    # Index the first anchor by ID
    anchor_rows = dict()

    a1_by_id = dict()
    for sa_1 in a1.strands:
        a1_by_id[sa_1.idstr] = sa_1
    
    for sa_2 in a2.strands:
        if sa_2.idstr in a1_by_id:
            
            # This read strand is anchored on both ends
            sa_1 = a1_by_id[sa_2.idstr]

            out.append( (sa_1, sa_2) )

    return out
            
def add_read_state_from_anchor_strands(lib, sa_1, sa_2):
    
    start_ei = sa_1.event_idx
    end_ei = sa_2.event_idx
    stride = 1

    if start_ei > end_ei:
        stride = -1

    strand_idx = 0
    if sa_1.strand == 'c':
        strand_idx = 1

    rs = CReadStateInterface(sa_1.sr_idx, start_ei, end_ei, strand_idx, stride, sa_1.rc)
    lib.add_read_state(rs)

# Import squiggle reads into C code
def load_reads(lib, squiggle_reads):

    for sr in squiggle_reads:
        if sr is None:
            continue
        
        pm_params = []
        event_params = []

        # Setup pore model and event sequence data for each strand of the read
        for s in ('t', 'c'):

            # Pore Model
            level_mean = sr.pm[s].model_level_mean.ctypes.data_as(c_double_p)
            level_stdv = sr.pm[s].model_level_stdv.ctypes.data_as(c_double_p)
            sd_mean = sr.pm[s].model_sd_mean.ctypes.data_as(c_double_p)
            sd_stdv = sr.pm[s].model_sd_stdv.ctypes.data_as(c_double_p)
            
            scale = sr.pm[s].scale
            shift = sr.pm[s].shift
            drift = sr.pm[s].drift
            var = sr.pm[s].var
            
            pm_params.append(CPoreModelInterface(1024, level_mean, level_stdv, sd_mean, sd_stdv, scale, shift, drift, var))

            # Events
            n_events = len(sr.event_level[s])
            level = sr.event_level[s].ctypes.data_as(c_double_p)
            stdv = sr.event_stdv[s].ctypes.data_as(c_double_p)
            time = sr.event_time[s].ctypes.data_as(c_double_p)
            event_params.append(CEventSequenceInterface(n_events, level, stdv, time))

        #
        params = CSquiggleReadInterface((pm_params[0], pm_params[1]), 
                                        (event_params[0], event_params[1]))
        lib.add_read(params)



#
# Build the map from read indices to fast5 files
#
fast5_fofn_fn = 'r73.map.fofn'
fast5_fofn_fh = open(fast5_fofn_fn)

f5_files = []
for l in fast5_fofn_fh:
    f5_files.append(l.rstrip())

#
# Load squiggle reads
#
use_poabase_signals = False
do_training = False

clustal_filename = "clustal-0.out"
clustal = Clustal(clustal_filename)
#cons_row = clustal.get_consensus_row()
print "!!!!!!!! using first read as consensus row !!!!!!!"

cons_row = 0
read_rows = clustal.get_read_rows()

# Initialize traversal
n_reads = len(read_rows)

# Parse the clustal row names to generate POAReads
poa_reads = list()
row_to_poa_read_idx = dict()
for rr in read_rows:
    pa = unpack_poa_id(clustal.alignment[rr].id)

    # Append read and update map
    poa_reads.append(pa)
    
    idx = len(poa_reads) - 1
    row_to_poa_read_idx[rr] = idx

sys.stderr.write("Loading squigglereads\n")
squiggle_reads = list()
poa_idx_to_sr_idx = dict()

seen_ids = dict()

for (poa_idx, pr) in enumerate(poa_reads):

    # A read may appear in multiple rows of the MA
    # if the aligner emitted a few local alignments
    # Only load a squiggleread once
    if pr.read_id in seen_ids:
        poa_idx_to_sr_idx[poa_idx] = seen_ids[pr.read_id]
        continue

    print 'Loading', pr.read_id, f5_files[pr.read_id]
    squiggle_reads.append(SquiggleRead(f5_files[pr.read_id]))

    sr_idx = len(squiggle_reads) - 1
    poa_idx_to_sr_idx[poa_idx] = sr_idx
    seen_ids[pr.read_id] = sr_idx

sys.stderr.write("Done loading squigglereads\n")


#
# Build anchors
#

DEBUG=False
ANCHOR_DISTANCE = 50 # in bases along the consensus row

# this list stores the number of actual bases seen in each row up to current column
n_bases = [0] * n_reads
anchors = list()

cons_kidx = 0
consensus_sequence = str(clustal.alignment[cons_row].seq)

n_cols = len(clustal.alignment[cons_row])
col = 0

#
while col < n_cols:

    # Is there a consensus base in this column?
    cons_base = clustal.alignment[cons_row][col]
    cons_kmer = clustal.get_kmer(cons_row, col, 5)

    # Should we make an anchor here?
    if cons_kidx % ANCHOR_DISTANCE == 0 and cons_base != '-' and len(cons_kmer) == 5:

        print 'ANCHOR', cons_kidx
        current_anchor = Anchor(col, cons_kidx)

        # Make anchor
        for rr in read_rows:
            b = clustal.alignment[rr][col]

            # Read metadata
            poa_read = poa_reads[rr]

            if poa_read.is_base and not use_poabase_signals:
                continue

            poa_idx = row_to_poa_read_idx[rr]
            sr_idx = poa_idx_to_sr_idx[poa_idx]
            sr = squiggle_reads[sr_idx]

            read_kidx = n_bases[rr] + poa_read.start
            
            last_read_kidx = poa_read.stop - 5
            if poa_read.is_base:
                last_read_kidx = sr.get_2D_length() - 5

            print "\tROW", rr, poa_reads[rr].read_id, read_kidx, last_read_kidx

            # skip reads that are not aligned in this range
            if (read_kidx <= poa_read.start and b == '-') or (read_kidx >= last_read_kidx and b == '-'):
                continue

            skmer = cons_kmer

            t_rc = 0
            c_rc = 1
            if poa_read.strand == 'c':

                # switch kmer index and kmer sequence to opposite strand
                read_kidx = sr.flip_k_idx_strand(read_kidx, 5)
                skmer = revcomp(skmer)

                t_rc = 1
                c_rc = 0

            # calculate closest events to this position
            (tdiff, ti) = calculate_diff_to_expected(sr, read_kidx, 't', skmer)
            (cdiff, ci) = calculate_diff_to_expected(sr, read_kidx, 'c', revcomp(skmer))
            
            # Fail if both events couldn't be found
            if ti == -1 or ci == -1:
                continue

            for (ei, strand, rc, diff) in ( (ti, 't', t_rc, tdiff), (ci, 'c', c_rc, cdiff)):
                strand_anchor = StrandAnchor(str(sr_idx) + "/" + strand, rr, sr_idx, strand, ei, rc, diff)
                current_anchor.add_strand(strand_anchor)
                current_anchor.n_samples += 1
                current_anchor.sum_diff += abs(diff)

        anchors.append(current_anchor)

    #
    # Update indices
    #
    if cons_base != '-':
        cons_kidx += 1

    # Update base counts
    for rr in read_rows:
        b = clustal.alignment[rr][col]
        if b != '-':
            n_bases[rr] += 1
    col += 1

print "!!!!!!! move last anchor to end of consensus !!!!!!"

#
# Initialize library
#
lib_hmmcons_fast = cdll.LoadLibrary("../nanopolish/hmmcons_fast.so")
lib_hmmcons_fast.initialize()

#
# Add reads
#
load_reads(lib_hmmcons_fast, squiggle_reads)

#
# Call consensus
#

# Add anchors to the library
for i in xrange(0, len(anchors)):
    
    a = anchors[i]
    
    anchor_rows = dict()

    # Index anchors by squiggle read idx and strand
    anchors_by_id = dict()
    for sa in a.strands:
        anchors_by_id[sa.idstr] = sa

    # Tell the library we are starting to add anchors
    lib_hmmcons_fast.start_anchored_column()

    # We pass 2 read anchors to the library for every read
    # event if it doesn't have events here. This must be done
    # in the same order that we passed reads to the library
    for (sr_idx, sr) in enumerate(squiggle_reads):
        for strand in ('t', 'c'):
            idstr = str(sr_idx) + "/" + strand

            ra = CReadAnchorInterface(-1, 0)
            if idstr in anchors_by_id:
                sa = anchors_by_id[idstr]
                ra = CReadAnchorInterface(sa.event_idx, sa.rc)
                anchor_rows[sa.row_id] = 1

            lib_hmmcons_fast.add_read_anchor(ra)

    # Add sequences to this segment if its not the last
    if i != len(anchors) - 1:
        next_anchor = anchors[i + 1]
        base_sequence = clustal.get_sequence_plus_k(cons_row, a.column, next_anchor.column, 5)
        lib_hmmcons_fast.add_base_sequence(base_sequence)
    
        for row in anchor_rows:
            read_sub = clustal.get_sequence_plus_k(row, a.column, next_anchor.column, 5)
            lib_hmmcons_fast.add_alt_sequence(read_sub)

    # Tell the library we are finished adding data for this segment
    lib_hmmcons_fast.end_anchored_column()

lib_hmmcons_fast.run_splice()

sys.exit(1)

#
# Step 1. Learn parameters of the model
#
if do_training:
    for i in xrange(0, len(anchors) - 1):
        lib_hmmcons_fast.clear_state()

        a1 = anchors[i]
        a2 = anchors[i+1]

        # Add the consensus sequence as the initial candidate consensus
        candidate_sub = clustal.get_sequence_plus_k(cons_row, a1.column, a2.column, 5)
        
        lib_hmmcons_fast.add_candidate_consensus(candidate_sub)
        
        # Initialize the c library using the reads that are anchored here
        anchored_strand_pairs = pair_anchor_strands(a1, a2)

        for (sa_1, sa_2) in anchored_strand_pairs:
            add_read_state_from_anchor_strands(lib_hmmcons_fast, sa_1, sa_2)
        
        lib_hmmcons_fast.learn_segment()

    lib_hmmcons_fast.train()

# Call a consensus between the first two anchors
original_consensus = ""
fixed_consensus = ""

for i in xrange(0, len(anchors) - 1):
    lib_hmmcons_fast.clear_state()

    a1 = anchors[i]
    a2 = anchors[i+1]

    if DEBUG:
        a1 = anchors[3]
        a2 = anchors[4]

    # Add the consensus sequence as the initial candidate consensus
    candidate_sub = clustal.get_sequence_plus_k(cons_row, a1.column, a2.column, 5)

    if DEBUG and False:
        truth = "ATTGCGCGATCAGTTGCCTGCCACCACTGTTCCGGGTCTTGTTCCGACCAGAGTGGATGCGGGCGCGAAACGGTCAGCTTTTCCGTTTGCGCAGC"
        lib_hmmcons_fast.add_candidate_consensus(truth)
    else:
        lib_hmmcons_fast.add_candidate_consensus(candidate_sub)
    
    print "CONSENSUS -- %d to %d %s" % (a1.column, a2.column, candidate_sub)

    # Set up event sequences
    anchor_rows = dict()

    anchored_strand_pairs = pair_anchor_strands(a1, a2)
    
    for (sa_1, sa_2) in anchored_strand_pairs:
        add_read_state_from_anchor_strands(lib_hmmcons_fast, sa_1, sa_2)
        anchor_rows[sa_1.row_id] = 1

    # Add the sequences of each row as alternates
    for row in anchor_rows:
        read_sub = clustal.get_sequence_plus_k(row, a1.column, a2.column, 5)
        lib_hmmcons_fast.add_candidate_consensus(read_sub)

    lib_hmmcons_fast.run_splice()

    # extract results
    lib_hmmcons_fast.get_consensus_result.restype = c_char_p
    result = lib_hmmcons_fast.get_consensus_result()
    print "POACON[%d]: %s" % (i, candidate_sub)
    print "RESULT[%d]: %s" % (i, result)
    
    if original_consensus == "":
        original_consensus = candidate_sub
        fixed_consensus = result
    else:
        original_consensus = original_consensus + candidate_sub[5:]
        fixed_consensus = fixed_consensus + result[5:]
    
    print "ORIGINAL: ", original_consensus
    print "FIXED: ", fixed_consensus

    if DEBUG:
        sys.exit(0)

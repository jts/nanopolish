from SquiggleRead import *
from NanoUtils import *
from Clustal import *
from ctypes import *
from collections import defaultdict, namedtuple

ReadAnchor = namedtuple('ReadAnchor', ['row_id', 'base'])

class Anchor:
    def __init__(self, column, cons_base):
        self.column = column
        self.cons_base = cons_base
        self.reads = list()

    def add_read(self, ra):
        self.reads.append(ra)

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

# Get the closest event to the indicated
# 2D kmer for the given squiggle read
def get_closest_event_to(sr, kidx, strand):
    kmer = sr.get_2D_kmer_at(kidx, 5)
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

    # evaluate levels
    b_level = sr.get_drift_corrected_event_level(first_event_before, strand)
    a_level = sr.get_drift_corrected_event_level(first_event_after, strand)
    #print 'SEARCH', first_event_before, first_event_after, b_level, a_level, e_level

    if abs(b_level - e_level) < abs(a_level - e_level):
        return first_event_before
    else:
        return first_event_after

#
# 
#
lib_hmmcons_fast = cdll.LoadLibrary("../nanopolish/hmmcons_fast.so")
lib_hmmcons_fast.initialize()

# Convenience type for working with ctypes
c_double_p = POINTER(c_double)

class CPoreModelInterface(Structure):
    _fields_ = [('n_states', c_int),
                ('pore_model_mean', c_double_p),
                ('pore_model_sd', c_double_p),
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

#
# Load reads
# 
test_data = [ (0,    39,   'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch101_file106_strand.fast5", 32),
              (1608, 1650, 'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch289_file43_strand.fast5", 5106),
              (386,  425,  'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch106_file136_strand.fast5", 198),
              (2747, 2784, 'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch364_file36_strand.fast5", 6852) ]

#
# Old test code
#

#reads = []
#for (start, stop, strand, fn) in test_data:
#    reads.append(SquiggleRead(fn))

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
clustal_filename = "clustal-32.out"
ca = Clustal(clustal_filename)
cons_row = ca.get_consensus_row()
read_rows = ca.get_read_rows()

# Initialize traversal
n_reads = len(read_rows)

# Parse the clustal row names to generate POAReads
poa_reads = list()
row_to_poa_read_idx = dict()
for rr in read_rows:
    pa = unpack_poa_id(ca.alignment[rr].id)

    # Append read and update map
    poa_reads.append(pa)
    
    idx = len(poa_reads) - 1
    row_to_poa_read_idx[rr] = idx

sys.stderr.write("Loading squigglereads\n")
squiggle_reads = list()
poa_idx_to_sr_idx = dict()

seen_ids = dict()

for (poa_idx, pr) in enumerate(poa_reads):

    #
    # DEBUG: only load a subset
    #
    in_debug_set = False
    for t in test_data:
        if pr.read_id == t[4]:
            in_debug_set = True

    #if not in_debug_set:
    #    continue

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
# Initialize squiggle reads in C code
#
for sr in squiggle_reads:
    if sr is None:
        continue
    
    pm_params = []
    event_params = []

    # Setup pore model and event sequence data for each strand of the read
    for s in ('t', 'c'):

        # Pore Model
        mean = sr.pm[s].model_mean.ctypes.data_as(c_double_p)
        sd = sr.pm[s].model_sd.ctypes.data_as(c_double_p)
        scale = sr.pm[s].scale
        shift = sr.pm[s].shift
        drift = sr.pm[s].drift
        var = sr.pm[s].var
        pm_params.append(CPoreModelInterface(1024, mean, sd, scale, shift, drift, var))

        # Events
        n_events = len(sr.event_level[s])
        level = sr.event_level[s].ctypes.data_as(c_double_p)
        stdv = sr.event_stdv[s].ctypes.data_as(c_double_p)
        time = sr.event_time[s].ctypes.data_as(c_double_p)
        event_params.append(CEventSequenceInterface(n_events, level, stdv, time))

    #
    params = CSquiggleReadInterface((pm_params[0], pm_params[1]), 
                                    (event_params[0], event_params[1]))
    lib_hmmcons_fast.add_read(params)

# Build anchors

# compute anchor positions
cons_kidx = 0

# this list stores the number of actual bases seen in each row
# up to the current column
n_bases = [0] * n_reads
anchors = list()

distance_to_last_anchor = 0
MIN_DISTANCE = 50
MAX_DISTANCE = 100

consensus_sequence = str(ca.alignment[cons_row].seq)

for col in xrange(0, 2000):
    cons_base = ca.alignment[cons_row][col]

    if cons_base != '-' and (len(anchors) == 0 or distance_to_last_anchor >= MIN_DISTANCE):

        # Should we anchor here?
        use_anchor = False

        #cons_kmer = consensus_sequence[cons_kidx:cons_kidx + 5]
        cons_kmer = ca.get_kmer(cons_row, col, 5)
        print 'ANCHOR', cons_kidx, cons_kmer

        max_res = 0
        n_evidence = 0
        for rr in read_rows:
            b = ca.alignment[rr][col]

            poa_read = poa_reads[rr]
            read_kidx = n_bases[rr] + poa_reads[rr].start

            # skip reads that are not aligned in this range
            if not poa_read.is_base and (read_kidx <= poa_read.start or read_kidx >= poa_read.stop):
                continue

            read_kidx = n_bases[rr] + poa_reads[rr].start
            poa_idx = row_to_poa_read_idx[rr]
            sr_idx = poa_idx_to_sr_idx[poa_idx]
            sr = squiggle_reads[sr_idx]
            orientation = poa_reads[rr].strand

            print "\tROW", rr, poa_idx, sr_idx, poa_reads[rr].read_id

            do_template_rc = 0
            do_complement_rc = 1

            if orientation == 'c':
                read_kidx = sr.flip_k_idx_strand(read_kidx, 5)
                do_template_rc = 1
                do_complement_rc = 0

            ti = get_closest_event_to(sr, read_kidx, 't')
            ci = get_closest_event_to(sr, read_kidx, 'c')

            #print "EVENTS", poa_reads[rr].read_id, ti, ci
            for (ei, strand, do_rc) in [ (ti, 't', do_template_rc), (ci, 'c', do_complement_rc) ]:
            
                if ei == -1:
                    continue

                # Observed event   
                level = sr.get_drift_corrected_event_level(ei, strand)

                # Expected event
                kmer = cons_kmer
                if do_rc:
                    kmer = revcomp(cons_kmer)

                k_level = sr.get_expected_level(kmer, strand)
                k_sd = sr.get_expected_sd(kmer, strand)
                res = (level - k_level) / k_sd

                print '\t\tRESIDUAL', strand, level, kmer, k_level, res
                if abs(res) > abs(max_res):
                    max_res = res
                n_evidence += 1
        if abs(max_res) < 3 and n_evidence >= 4:
            use_anchor = True
                    
        # Build an anchor here
        if use_anchor:
            print "BUILD", cons_kidx, cons_base
            
            anchor = Anchor(col, cons_kidx)

            for rr in read_rows:
                read_kidx = n_bases[rr] + poa_reads[rr].start
            
                poa_read = poa_reads[rr]
                read_kidx = n_bases[rr] + poa_reads[rr].start

                if poa_read.is_base or (read_kidx > poa_read.start and read_kidx < poa_read.stop):
                    read_anchor = ReadAnchor(rr, read_kidx)
                    anchor.add_read(read_anchor)
    
            anchors.append(anchor)
            distance_to_last_anchor = 0

    distance_to_last_anchor += 1
    #
    # Update indices
    #
    if cons_base != '-':
        cons_kidx += 1

    # Update base counts
    for rr in read_rows:
        b = ca.alignment[rr][col]
        if b != '-':
            n_bases[rr] += 1

# Call a consensus between the first two anchors
original_consensus = ""
fixed_consensus = ""

for i in xrange(len(anchors) - 1):
    lib_hmmcons_fast.clear_state()

    a1 = anchors[i]
    a2 = anchors[i+1]

    # Add the consensus sequence as the initial candidate consensus
    consensus = ca.alignment[cons_row].seq
    candidate_sub = str(consensus[a1.column:a2.column+5]).replace('-', '')
    lib_hmmcons_fast.add_candidate_consensus(candidate_sub)

    # Index the first anchor by read row
    a1_by_row = dict()
    for ra_1 in a1.reads:
        a1_by_row[ra_1.row_id] = ra_1

    for ra_2 in a2.reads:
        if ra_2.row_id in a1_by_row:
            
            # This row is anchored on both ends
            ra_1 = a1_by_row[ra_2.row_id]

            # Build a ReadState object to pass to C
            k_start_idx = ra_1.base
            k_stop_idx = ra_2.base
            read_sub = str(ca.alignment[ra_2.row_id].seq[a1.column:a2.column+5]).replace('-', '')
            lib_hmmcons_fast.add_candidate_consensus(read_sub)
            
            orientation = poa_reads[ra_1.row_id].strand
            poa_idx = row_to_poa_read_idx[ra_1.row_id]
            sr_idx = poa_idx_to_sr_idx[poa_idx]

            print 'READ_INDEX', ra_1.row_id, poa_idx, sr_idx

            sr = squiggle_reads[sr_idx]

            if orientation == 'n':
                t_stride = 1
                c_stride = -1
                t_rc = 0
                c_rc = 1
            else:
                t_stride = -1
                c_stride = 1
                t_rc = 1
                c_rc = 0
                k_start_idx = sr.flip_k_idx_strand(k_start_idx, 5)
                k_stop_idx = sr.flip_k_idx_strand(k_stop_idx, 5)

            k1 = sr.get_2D_kmer_at(k_start_idx, 5)
            k2 = sr.get_2D_kmer_at(k_stop_idx, 5)

            print '\tEVENTS', orientation, k_start_idx, k_stop_idx, k1, k2, 'START', sr.event_map['2D'][k_start_idx], 'END', sr.event_map['2D'][k_stop_idx]

            if orientation == 'c':
                k1 = revcomp(k1)
                k2 = revcomp(k2)

            t_start_ei = get_closest_event_to(sr, k_start_idx, 't')
            c_start_ei = get_closest_event_to(sr, k_start_idx, 'c')
            t_stop_ei = get_closest_event_to(sr, k_stop_idx, 't')
            c_stop_ei = get_closest_event_to(sr, k_stop_idx, 'c')

            if t_start_ei != -1 and t_stop_ei != -1:
                print '\tT_EVENTS', t_start_ei, t_stop_ei, t_stride
                t_rs = CReadStateInterface(sr_idx, t_start_ei, t_stop_ei, 0, t_stride, t_rc)
                lib_hmmcons_fast.add_read_state(t_rs)
            
            if c_start_ei != -1 and c_stop_ei != -1:
                print '\tC_EVENTS', c_start_ei, c_stop_ei, c_stride
                c_rs = CReadStateInterface(sr_idx, c_start_ei, c_stop_ei, 1, c_stride, c_rc)
                lib_hmmcons_fast.add_read_state(c_rs)

    #lib_hmmcons_fast.run_selection()
    #lib_hmmcons_fast.run_mutation()
    lib_hmmcons_fast.run_rewrite()
    #lib_hmmcons_fast.run_consensus()

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
sys.exit(0)

#
# Initialize HMM by telling C code where the events
# start and end for each read
#
for (ri, sr) in enumerate(reads):

    k_start_idx = test_data[ri][0]
    k_stop_idx = test_data[ri][1]
    orientation = test_data[ri][2]

    if orientation == 'n':
        t_stride = 1
        c_stride = -1
        t_rc = 0
        c_rc = 1
    else:
        t_stride = -1
        c_stride = 1
        t_rc = 1
        c_rc = 0
        k_start_idx = reads[ri].flip_k_idx_strand(k_start_idx, 5)
        k_stop_idx = reads[ri].flip_k_idx_strand(k_stop_idx, 5)

    k1 = reads[ri].get_2D_kmer_at(k_start_idx, 5)
    k2 = reads[ri].get_2D_kmer_at(k_stop_idx, 5)

    if orientation == 'c':
        k1 = revcomp(k1)
        k2 = revcomp(k2)

    (t_start_ei, c_start_ei) = reads[ri].event_map['2D'][k_start_idx][0]
    (t_stop_ei, c_stop_ei) = reads[ri].event_map['2D'][k_stop_idx][0]
    
    print k1, k2, t_start_ei, t_stop_ei, c_start_ei, c_stop_ei, t_stride, c_stride

    t_rs = CReadStateInterface(ri, t_start_ei, t_stop_ei, 0, t_stride, t_rc)
    lib_hmmcons_fast.add_read_state(t_rs)
    
    c_rs = CReadStateInterface(ri, c_start_ei, c_stop_ei, 1, c_stride, c_rc)
    lib_hmmcons_fast.add_read_state(c_rs)

lib_hmmcons_fast.run_debug()
#lib_hmmcons_fast.run_mutation()
#lib_hmmcons_fast.run_consensus()
sys.exit(0)

#
# OLD
#

test_data = [ (0,    'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch101_file106_strand.fast5"),
              (1608, 'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch289_file43_strand.fast5"),
              (386,  'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch106_file136_strand.fast5"),
              (2747, 'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch364_file36_strand.fast5") ]

reads = []
for (offset, strand, fn) in test_data:
    reads.append(SquiggleRead(fn))

kmers = [ "AACAG" ]
paths = generate_kmer_paths(kmers[-1])
rc_paths = []
for p in paths:
    rc_paths.append(revcomplist(p))

for (p, rc_p) in zip(paths, rc_paths):

    for (ri, r) in enumerate(reads):

        # skip for now
        if test_data[ri][1] == 'c':
            continue
    
    (t_ei, c_ei) = reads[ri].event_map['2D'][0][0]
    print t_ei, c_ei

    #    
    # Read 1 example
    #
    
    # need to handle opposite strand here, tricky
    if False:
        ri = 1
        k_idx = test_data[ri][0]
        k_idx = reads[ri].flip_k_idx_strand(k_idx, 5)
        print test_data[ri][0], k_idx
        (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]
        f_t_1 = reads[ri].hmm2(rc_p, 't', t_ei, -1)
        f_c_1 = reads[ri].hmm2(p, 'c', c_ei, 1)
        print 'FINAL', f_t_1, f_c_1, p
        continue

    #for s in xrange(-3, 3):
    #    kmer = reads[ri].get_2D_kmer_at(k_idx + s, 5)
    #    (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx + s][0]
    #    print s, kmer, revcomp(kmer), t_ei, c_ei, [(shortevent(x), shortevent(y)) for (x,y) in reads[ri].get_events_for_2D_kmer(k_idx + s)]
    #break

    #
    # Read 2 example
    #
    if False:
        ri = 2
        k_idx = test_data[ri][0]
        (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]
        f_t_2 = reads[ri].hmm2(p, 't', t_ei, 1)
        f_c_2 = reads[ri].hmm2(rc_p, 'c', c_ei, -1)
        #print f_t_2, f_c_2, p
        continue

    #
    # Read 0 example
    #
    if True:
        ri = 0
        k_idx = 0
        (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]

        # template events
        f_t_0 = reads[ri].hmm2(p, 't', t_ei, 1)
        f_c_0 = reads[ri].hmm2(rc_p, 'c', c_ei, -1)
        #print f_t_0 + f_c_0 + f_t_2 + f_c_2, f_t_0, f_c_0, f_t_2, f_c_2, p
        continue

    #
    base_index = test_data[2][0]
    events = reads[2].event_map['2D'][base_index]
    offset = events[0][0]

    f2 = 0
    if offset != -1:
        f2 = reads[2].hmm2(e, 't', offset)
    
    s = f1 + f2
    print s, f1, f2, e
    #print e, f2

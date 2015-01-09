from SquiggleRead import *
from NanoUtils import *
from ctypes import *
from collections import defaultdict

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
test_data = [ (0,    39,   'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch101_file106_strand.fast5"),
              (1608, 1650, 'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch289_file43_strand.fast5"),
              (386,  425,  'c', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch106_file136_strand.fast5"),
              (2747, 2784, 'n', "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch364_file36_strand.fast5") ]

reads = []
for (start, stop, strand, fn) in test_data:
    reads.append(SquiggleRead(fn))

#
# Pass reads into C code
#
for (ri, sr) in enumerate(reads):
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

transition_counts = defaultdict(int)

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

from SquiggleRead import *
from NanoUtils import *

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

    #    
    # Read 1 example
    #
    
    # need to handle opposite strand here, tricky
    ri = 1
    k_idx = test_data[ri][0]
    k_idx = reads[ri].flip_k_idx_strand(k_idx, 5)
    print test_data[ri][0], k_idx
    (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]
    f_t_1 = reads[ri].hmm2(rc_p, 't', t_ei, -1)
    f_c_1 = reads[ri].hmm2(p, 'c', c_ei, 1)
    print 'FINAL', f_t_1, f_c_1, p

    #for s in xrange(-3, 3):
    #    kmer = reads[ri].get_2D_kmer_at(k_idx + s, 5)
    #    (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx + s][0]
    #    print s, kmer, revcomp(kmer), t_ei, c_ei, [(shortevent(x), shortevent(y)) for (x,y) in reads[ri].get_events_for_2D_kmer(k_idx + s)]
    #break
    continue

    #
    # Read 2 example
    #
    ri = 2
    k_idx = test_data[ri][0]
    (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]
    f_t_2 = reads[ri].hmm2(p, 't', t_ei)
    f_c_2 = reads[ri].hmm2(rc_p, 'c', c_ei)
    #print f_t_2, f_c_2, p
    #continue

    #
    # Read 0 example
    #
    ri = 0
    k_idx = 0
    (t_ei, c_ei) = reads[ri].event_map['2D'][k_idx][0]

    # template events
    f_t_0 = reads[ri].hmm2(p, 't', t_ei)
    f_c_0 = reads[ri].hmm2(rc_p, 'c', c_ei)
    print f_t_0 + f_c_0 + f_t_2 + f_c_2, f_t_0, f_c_0, f_t_2, f_c_2, p
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

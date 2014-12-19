import itertools
from poretools import Fast5FileSet, Fast5File
from NanoUtils import *
from PoreModel import *
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import collections  as mc

# A SquiggleRead represents a nanopore read
# as a list of events and a mapping between
# template and complement events
class SquiggleRead:

    def __init__(self, filename):

        self.fh = Fast5File(filename)
        self.load_events()
        self.load_poremodel()

    def load_events(self):
        self.events = {}
        self.events['t'] = self.fh.get_template_events()
        self.events['c'] = self.fh.get_complement_events()

        self.event_map = {}
        self.event_map['t'] = self.build_1D_event_map('t')
        self.event_map['c'] = self.build_1D_event_map('c')
        self.event_map['2D'] = self.build_2D_event_map()

    def load_poremodel(self):
        self.pm = {}
        self.pm['t'] = PoreModel(self.fh, 'template')
        self.pm['c'] = PoreModel(self.fh, 'complement')

        #self.pm['t'].plot("template.png")
        #self.pm['c'].plot("complement.png")

    def has_events(self):
        return len(self.events['t']) + len(self.events['c']) > 0

    def get_template_sequence(self):
        return self.fh.get_fastas('fwd')[0].seq
    
    def get_complement_sequence(self):
        return self.fh.get_fastas('rev')[0].seq
    
    def get_2D_sequence(self):
        return self.fh.get_fastas('2D')[0].seq

    def get_2D_kmer_at(self, k_idx, k):
        return self.fh.get_fastas('2D')[0].seq[k_idx:k_idx+k]

    # Build a map from called read k-mer to template/complement event
    def build_2D_event_map(self):

        try:
            twod_alignment = list(self.fh.hdf5file['/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment'])
            
            seq = self.get_2D_sequence()
            event_map = list()
            k_idx = 0
            
            for (tidx, cidx, kmer) in twod_alignment:

                # Compute the index of this kmer in the 2D sequence
                while k_idx < len(seq) and kmer != seq[k_idx:k_idx+5]:
                    k_idx += 1

                # Append new lists to the list
                while len(event_map) <= k_idx:
                    event_map.append(list())

                event_map[k_idx].append((tidx, cidx))
            return event_map

        except Exception, e:
            print e
            return None
               
    # Map from a sequence k-mer to the index of the events showing it
    def build_1D_event_map(self, strand):
        event_map = list()

        events = self.events[strand]
        e_idx = 0
        k_idx = 0

        for (e_idx, e) in enumerate(events):
            k_idx += e.move
            while len(event_map) <= k_idx:
                event_map.append(list())

            event_map[k_idx].append(e_idx)
        return event_map

    #
    def get_events_for_kmer(self, k_idx, strand):

        out = []
        for e_idx in self.event_map[strand][k_idx]:
            out.append(self.events[strand][e_idx])
        return out

    # Get the events for the 2D basecall for k-mer k_idx
    def get_events_for_2D_kmer(self, k_idx):

        out = []

        for (tidx, cidx) in self.event_map['2D'][k_idx]:
            te = None
            ce = None
            if tidx >= 0:
                te = self.events['t'][tidx]
            if cidx >= 0:
                ce = self.events['c'][cidx]
            
            out.append((te, ce))
        return out
    
    # Get the events immediately before the 2D basecall at k-mer k_idx
    def get_event_before_2D_kmer(self, k_idx):

        out = []

        (tidx, cidx) = self.event_map['2D'][k_idx][0]

        te = None
        ce = None

        if tidx >= 0:
           te = self.events['t'][tidx - 1]
        if cidx >= 0:
           ce = self.events['c'][cidx + 1]
        return (te, ce)
    
    # Get the events immediately following the 2D basecall at k-mer k_idx
    def get_event_after_2D_kmer(self, k_idx):

        out = []
        (tidx, cidx) = self.event_map['2D'][k_idx][-1]
        te = None
        ce = None
        if tidx >= 0:
            te = self.events['t'][tidx + 1]
        if cidx >= 0:
            ce = self.events['c'][cidx - 1]
        return (te,ce)

    def get_expected_level(self, k_mer, strand):
        return self.pm[strand].get_expected(k_mer)

    def format_kmer_event_pair(self, te, ce, expected_kmer):
        tstr = self.format_kmer_event(te, expected_kmer, 't')
        cstr = self.format_kmer_event(ce, revcomp(expected_kmer), 'c')
        return expected_kmer + "\n" + tstr + cstr

    def format_kmer_event(self, e, expected_kmer, strand):
        ostr = "No event\n"
        if e != None:
            el = self.get_expected_level(expected_kmer, strand)
            d = e.mean - el 
            ostr = "Level: %.1f Called: %s Expected: %.1f D: %.1f\n" % (e.mean, e.model_state, el, d)
        return ostr

    # Calculate the coordinate of this k-mer on the other strand
    def flip_k_idx_strand(self, k_idx, k):
        return len(self.fh.get_fastas('2D')[0].seq) - k_idx - k

    def create_trace(self, m, n):
        out = list()
        for i in xrange(0, m):
            out.append(list(('-') * n))
        return out

    def hmm_align(self, seq):

        # Clamp lengths to 100!
        n_rows = 9 # number of events
        n_cols = 11 # number of k-mers

        strand = 't'

        # Break the sequence into k-mers for easy access
        kmers = []
        for i in xrange(0, len(seq)):
            kmers.append(seq[i:i+5])

        # Easy access to the event stream
        events = self.events[strand]

        M = np.zeros(shape=(n_rows, n_cols))
        E = np.zeros(shape=(n_rows, n_cols))
        K = np.zeros(shape=(n_rows, n_cols))
        
        trace_M = self.create_trace(n_rows, n_cols)
        trace_E = self.create_trace(n_rows, n_cols)
        trace_K = self.create_trace(n_rows, n_cols)

        # Initialize
        M[0,0] = log(1.0)
        E[0,0] = log(0.0)
        K[0,0] = log(0.0)

        for row in xrange(1, n_rows):
            M[row, 0] = log(0.0)
            E[row, 0] = log(0.0)
            K[row, 0] = log(0.0)
            trace_E[row][0] = 'E'
            trace_M[row][0] = 'E'

        for col in xrange(1, n_cols):
            M[0, col] = log(0.0)
            E[0, col] = log(0.0)
            K[0, col] = log(0.0)
            trace_K[0][col] = 'K'
            trace_M[0][col] = 'K'
        
        print M
        print E
        print K

        # Pseudocounts to calculate normalized transition probability
        t_mat = np.matrix([[90., 5., 5.], [90., 10., 5.], [90, 5., 10.]])

        t_m_sum = sum(t_mat[0, :])
        t_e_sum = sum(t_mat[1, :])
        t_k_sum = sum(t_mat[2, :])

        print t_mat[0,0], t_m_sum, 

        # Fixed transition probabilities
        t_mm = log( t_mat[0,0] / t_m_sum )
        t_me = log( t_mat[0,1] / t_m_sum )
        t_mk = log( t_mat[0,2] / t_m_sum )
        
        t_em = log( t_mat[1,0] / t_e_sum )
        t_ee = log( t_mat[1,1] / t_e_sum )
        t_ek = log( t_mat[1,2] / t_e_sum )
        
        t_km = log( t_mat[2,0] / t_k_sum )
        t_ke = log( t_mat[2,1] / t_k_sum )
        t_kk = log( t_mat[2,2] / t_k_sum )
        
        print 'Match transitions:', exp(t_mm), exp(t_me), exp(t_mk)
        print 'EIns  transitions:', exp(t_em), exp(t_ee), exp(t_ek)
        print 'KIns  transitions:', exp(t_km), exp(t_ke), exp(t_kk)

        for col in xrange(1, n_cols):
            for row in xrange(1, n_rows):

                i = row - 1
                j = col - 1

                # Probability of matching k_i to e_j
                # This is P(k_i, e_j | M) = P(e_j | k_i, M) P(k_i | M)
                if i >= 0 and j >= 0:
                    l_p_m = self.pm[strand].get_log_probability(events[i].mean, kmers[j])
                else:
                    l_p_m = log(0.0)

                # Probability of an inserted event
                # We model this by calculating the probability that e_j
                # is actually the signal for the previous k-mer, k_(i-1)
                l_p_e = self.pm[strand].get_log_probability(events[i].mean, kmers[j])
                #l_p_e = log(0.1)

                # Probability of an inserted k-mer
                # THIS IS HACKY AND WRONG
                l_p_k = log(0.1)

                # Calculate M[i, j]
                d_m = t_mm + M[row - 1, col - 1]
                d_e = t_em + E[row - 1, col - 1]
                d_k = t_km + K[row - 1, col - 1]
                max_diag = max(d_m, d_e, d_k)

                # Fill trace matrix
                if max_diag == d_m:
                    trace_M[row][col] = 'M'
                elif max_diag == d_e:
                    trace_M[row][col] = 'E'
                elif max_diag == d_k:
                    trace_M[row][col] = 'K'

                M[row, col] = l_p_m + max_diag

                print 'M[%d, %d] = %.1f -- k_i: %s e_j: %.1f p_m: %.2f [%.1f %.1f %.1f] %c' % (row, col, M[row,col], kmers[j], events[i].mean, l_p_m, d_m, d_e, d_k, trace_M[row][col])

                # Calculate E[i, j]
                u_m = t_me + M[row - 1, col]
                u_e = t_ee + E[row - 1, col]
                u_k = t_ke + K[row - 1, col]

                max_up = max(u_m, u_e, u_k)

                if max_up == u_m:
                    trace_E[row][col] = 'M'
                elif max_up == u_e:
                    trace_E[row][col] = 'E'
                elif max_up == u_k:
                    trace_E[row][col] = 'K'

                E[row, col] = l_p_e + max_up

                print 'E[%d, %d] = %.1f -- k_i: %s e_j: %.1f p_x: %.2f [%.1f %.1f %.1f] %c' % (row, col, E[row,col], kmers[j], events[i].mean, l_p_e, u_m, u_e, u_k, trace_E[row][col])

                # Calculate K[i, j]
                l_m = t_mk + M[row, col - 1]
                l_k = t_kk + K[row, col - 1]
                l_e = t_ek + E[row, col - 1]
                max_left = max(l_m, l_k, l_e)

                if max_left == l_m:
                    trace_K[row][col] = 'M'
                elif max_left == l_k:
                    trace_K[row][col] = 'K'
                elif max_left == l_e:
                    trace_K[row][col] = 'E'

                K[row,col] = l_p_k + max_left
                print 'K[%d, %d] = %.1f -- k_i: %s e_j: %.1f p_y: %.2f [%.1f %1.f %.1f] %c' % (row, col, K[row,col], kmers[j], events[j].mean, l_p_k, l_m, l_e, l_k, trace_K[row][col])

        # Reconstruct path
        i = n_rows - 1
        j = n_cols - 1

        out = []
        curr_m = 'M'
        while i > 0 and j > 0:

            print 'Current matrix', curr_m
            print 'Current cell', i, j

            # What matrix was used to arrive at this cell of the current matrix?
            if curr_m == 'M':
                prev_m = trace_M[i][j]
            if curr_m == 'K':
                prev_m = trace_K[i][j]
            if curr_m == 'E':
                prev_m = trace_E[i][j]

            # Move to the previous cell based on the current matrix
            if curr_m == 'M':
                i -= 1
                j -= 1
                out.append('M')

            if curr_m == 'K':
                j -= 1
                out.append('K')

            if curr_m == 'E':
                i -= 1
                out.append('E')

            curr_m = prev_m    
        print ''.join(reversed(out))

# 
def example_hmm():
    # Textbook example
        
    # Initialize the model object
    model = Model( name="Rainy-Sunny" )

    # Initialize the two hidden states, with an appropriate discrete distribution
    rainy = State( DiscreteDistribution({ 'walk': 0.1, 'shop': 0.4, 'clean': 0.5 }), name='Rainy' )
    sunny = State( DiscreteDistribution({ 'walk': 0.6, 'shop': 0.3, 'clean': 0.1 }), name='Sunny' )

    # Add the states to the model
    model.add_state( rainy )
    model.add_state( sunny )

    # Now add the two transitions from the start of the model to the hidden states
    model.add_transition( model.start, rainy, 0.6 )
    model.add_transition( model.start, sunny, 0.4 )

    # Add the transitions from the hidden states to each other
    model.add_transition( rainy, rainy, 0.7 )
    model.add_transition( rainy, sunny, 0.3 )
    model.add_transition( sunny, rainy, 0.4 )
    model.add_transition( sunny, sunny, 0.6 )

    # Finalize the model structure
    model.bake( verbose=True )

    # observations
    day1 = [ 'walk', 'walk', 'shop', 'walk', 'clean', 'walk' ]

    print model.log_probability( day1 )

    print model.viterbi( day1 )
    print model.sample(10)

if __name__ == '__main__':

    test_file = "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch101_file106_strand.fast5"
    #test_file = "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch364_file36_strand.fast5"
    sr = SquiggleRead(test_file)
    
    #con0 = "AACAGTCCACTATTGGATGGTAAAGCCAACAGAAATTTTTACGCAAGCTAAAGCCCGGCAGATGATTATCTTGCCATATGACGTCAAACCGCGGTTTGAATGAAACGCTGGATGATATTTGCGAAGCATTGAGTATTATGT"
    #sr.hmm_align(con0)

    test_seq = sr.get_2D_sequence()[0:14]
    sr.hmm_align(test_seq)

    #sys.exit(1)

    for i in xrange(0, 10):
        kmer = sr.get_2D_kmer_at(i, 5)
        td_events = sr.get_events_for_2D_kmer(i)

        td_out = []
        for (x,y) in td_events:
            if x != None:
                x = (x.mean, x.model_state, sr.get_expected_level(kmer, 't'))
            td_out.append(x)

        print i, kmer, td_out


    if False:
        temp_seq = sr.get_template_sequence()
        comp_seq = sr.get_complement_sequence()
        twod_seq = sr.get_2D_sequence()
        
        temp_idx = 124
        t_events = sr.get_events_for_kmer(temp_idx, 't')

        # translate comp idx onto other strand
        comp_idx = 179
        comp_idx = len(comp_seq) - comp_idx - 5
        c_events = sr.get_events_for_kmer(comp_idx, 'c')
        
        twod_idx = 140
        twod_kmer = twod_seq[twod_idx:twod_idx+5]
        td_events = sr.get_events_for_2D_kmer(twod_idx)
        td_out = []
        for (x,y) in td_events:
            if x != None:
                x = (x.mean, x.model_state)
            if y != None:
                y = (y.mean, y.model_state, revcomp(y.model_state))
            td_out.append((x,y))
        
        print temp_idx, comp_idx, twod_idx

        print temp_seq[temp_idx:temp_idx+5], sr.get_expected_level(twod_kmer, 't'), [(x.mean, x.model_state, x.model_level) for x in t_events]
        print revcomp(comp_seq[comp_idx:comp_idx+5]), sr.get_expected_level(revcomp(twod_kmer), 'c'), [(x.mean, revcomp(x.model_state), x.model_level) for x in c_events]
        print twod_seq[twod_idx:twod_idx+5], td_out

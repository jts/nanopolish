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
            print tidx, cidx
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

if __name__ == '__main__':
    test_file = "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch101_file106_strand.fast5"
    #test_file = "../R73_data/downloads/LomanLabz_PC_Ecoli_K12_R7.3_2549_1_ch364_file36_strand.fast5"
    sr = SquiggleRead(test_file)

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

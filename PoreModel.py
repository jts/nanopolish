from SquiggleRead import *
from NanoUtils import *
from pylab import *

class PoreModel:
    def __init__(self, fh, strand):
        self.read_from_file(fh, strand)

    def read_from_file(self, fh, strand):
               
        obj = fh.hdf5file['/Analyses/Basecall_2D_000/BaseCalled_%s/Model' % (strand)]
        self.drift = obj.attrs['drift']
        self.scale = obj.attrs['scale']
        self.scale_sd = obj.attrs['scale_sd']
        self.shift = obj.attrs['shift']
        self.var = obj.attrs['var']
        self.var_2d = obj.attrs['var_sd']
        #print self.shift, self.drift, self.scale, self.var

        self.model = list(obj)
        self.rank_map = generate_mers_rank_map(dna_alphabet, 5)

    # Get the expected value for the given kmer
    def get_expected(self, kmer):
        # TODO: account for drift?
        return (self.model[self.rank_map[kmer]][2] + self.shift) * self.scale
    
    #
    def plot(self, filename):
        x = []
        y = []
        for (i, m) in enumerate(self.model):
            x.append(i)
            y.append(m[2])
        plt.plot(x, y)
        plt.savefig(filename, dpi=150)


from SquiggleRead import *
from NanoUtils import *
from pylab import *
from scipy.stats import norm

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
        self.var_sd = obj.attrs['var_sd']
        print self.shift, self.drift, self.scale, self.var

        self.model = np.array(obj)

        self.model_mean = np.array([x[2] for x in obj])
        self.model_sd = np.array([x[3] for x in obj])

        # coerece data into the right type to pass to the C HMM library
        self.model_mean = self.model_mean.astype(np.float64)
        self.model_sd = self.model_sd.astype(np.float64)

        self.rank_map = generate_mers_rank_map(dna_alphabet, 5)

    # Get the expected value for the given kmer
    def get_expected(self, kmer):
        # TODO: account for drift?
        return (self.model[self.rank_map[kmer]][2] + self.shift) * self.scale
    
    def get_sd(self, kmer):
        # TODO: account for drift?
        return self.model[self.rank_map[kmer]][3] * self.scale

    # Get the probability of seeing signal s, given the k-mer in the pore
    def get_probability(self, s, kmer):
        m = self.get_expected(kmer)
        sd = self.get_sd(kmer)
        p = norm.pdf(s, m, sd)
        print 'k: %s s: %.1f m: %.1f p: %.3f' % (kmer, m, s, p)
        return p
    
    # Get the probability of seeing signal s, given the k-mer in the pore
    def get_log_probability(self, s, kmer):
        m = self.get_expected(kmer)
        sd = self.get_sd(kmer)
        lp = norm.logpdf(s, m, sd)
        return lp
    
    #
    def plot(self, filename):
        x = []
        y = []
        for (i, m) in enumerate(self.model):
            x.append(i)
            y.append(m[2])
        plt.plot(x, y)
        plt.savefig(filename, dpi=150)


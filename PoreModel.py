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

        self.model = np.array(obj)

        self.model_level_mean = np.array([x['level_mean'] for x in obj])
        self.model_level_stdv = np.array([x['level_stdv'] for x in obj])
        
        self.model_sd_mean = np.array([x['sd_mean'] for x in obj])
        self.model_sd_stdv = np.array([x['sd_stdv'] for x in obj])

        # coerece data into the right type to pass to the C HMM library
        self.model_level_mean = self.model_level_mean.astype(np.float64)
        self.model_level_stdv = self.model_level_stdv.astype(np.float64)
        
        self.model_sd_mean = self.model_sd_mean.astype(np.float64)
        self.model_sd_stdv = self.model_sd_stdv.astype(np.float64)

        self.rank_map = generate_mers_rank_map(dna_alphabet, 5)

    def get_expected(self, kmer):
        # TODO: account for drift?
        #print kmer, self.rank_map[kmer], self.model[self.rank_map[kmer]][2], self.shift, self.scale
        return (self.model[self.rank_map[kmer]]['level_mean'] * self.scale) + self.shift
    
    def get_sd(self, kmer):
        # TODO: account for drift?
        return self.model[self.rank_map[kmer]]['level_stdv'] * self.var

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


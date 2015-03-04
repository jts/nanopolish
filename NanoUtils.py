import random
dna_alphabet = ['A', 'C', 'G', 'T']

# Generate a list of all k-mers
def generate_mers(alphabet, k):

    # initialze kmer list with an empty string
    kmers = list()
    kmers.append("")

    for i in range(0, k):

        # for every string in the output list,
        # generate a new string one base longer for
        # every symbol in the alphabet
        new_kmers = list()
        for o in kmers:
            for j in range(0, len(alphabet)):
                new_kmers.append(o + alphabet[j])
        kmers = new_kmers

    kmers.sort()
    return kmers

# Generate a map from k-mer -> lexicographic rank
def generate_mers_rank_map(alphabet, k):

    kmers = generate_mers(alphabet, k)

    # Build dictionary
    rank_dict = dict()
    for (i, o) in enumerate(kmers):
        rank_dict[o] = i
    return rank_dict

# Generate all of the kmers of the string
def str2kmers(s, k):
    out = []
    for i in xrange(0, len(s) - k + 1):
        out.append(s[i:i+k])
    return out

# reverse complement a sequence
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def revcomplist(l):
    out = []
    for k in l:
        out.append(revcomp(k))
    return out

def random_base():
    return dna_alphabet[random.randint(0, 3)]

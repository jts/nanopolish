#! /usr/bin/env python
# Generate a new model with a lower k-mer size than the input model

import sys
import argparse
import itertools
import numpy

alphabet = [ 'A', 'C', 'G','T' ]
def make_all_mers(k):
    return [ "".join(x) for x in itertools.product(alphabet, repeat=k) ]

parser = argparse.ArgumentParser( description='Reduce a 6-mer model to a 5-mer model')
parser.add_argument('-i', '--input', type=str, required=True)
args = parser.parse_args()

# Read the initial model from a file
f = open(args.input)

K = 0
model = dict()

header_lines_to_copy = { "#strand", "#kit", "#ont_model_name", "#alphabet" }
header_lines = list()
input_model_name = ""

for line in f:
    line = line.rstrip()
    fields = line.split()

    # copy then skip header lines
    if line[0] == '#' or line.find("kmer") == 0:
        if fields[0] in header_lines_to_copy:
            header_lines.append(line)
    else:
        # store the k-mer size
        if K == 0:
            K = len(fields[0])
        else:
            assert len(fields[0]) == K
        
        # store values 
        model[fields[0]] = tuple(fields[1:6])

# reduce the k-mer size by 1 and output the new model
P = K - 1
pmers = make_all_mers(P)

outname = args.input
kmer_str = str(K) + "mer"
pmer_str = str(P) + "mer"
assert(outname.find(kmer_str) != -1)
outname = args.input.replace(kmer_str, pmer_str)
outfile = open(outname, "w")
header_lines.append("#k\t" + str(P))
header_lines.append("#original_file\t" + args.input)

outfile.write("\n".join(header_lines) + "\n")
outfile.write("\t".join(["kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"]) + "\n")

num_samples_per_kmer = 1000

for pmer in pmers:
    
    kmers_with_pmer = [ pmer + a for a in alphabet ]
    samples = list()

    for kmer in kmers_with_pmer:
        
        # sample values from this gaussian
        samples += list(numpy.random.normal(model[kmer][0], model[kmer][1], num_samples_per_kmer))

    m = numpy.mean(samples)
    s = numpy.std(samples)
    out = [m, s, 0.0, 0.0, 0.0]
    outfile.write("\t".join([pmer] + [str(x) for x in out]) + "\n")

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
parser.add_argument('-t', '--type', type=str, required=False, default="dropmodel")
args = parser.parse_args()

# Read the initial model from a file
f = open(args.input)

K = 0
model = dict()

header_lines_to_copy = { "#strand", "#kit" }
header_lines = list()
input_model_name = ""

for line in f:
    line = line.rstrip()
    fields = line.split()

    # copy then skip header lines
    if line[0] == '#' or line.find("kmer") == 0:
        if fields[0] in header_lines_to_copy:
            header_lines.append(line)
        if fields[0] == "#model_name":
            input_model_name = fields[1]
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
print("#model_name\t" + args.input + ".dropmodel")
print("#type\t" + args.type)
print("\n".join(header_lines))
print("#derived_from\t" + input_model_name)

print "\t".join(["kmer", "level_mean", "level_stdv", "sd_mean", "sd_stdv"])

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
    print "\t".join([pmer] + [str(x) for x in out])

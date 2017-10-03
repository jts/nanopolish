#! /usr/bin/python

# This script takes a .model file provided by ONT and adds metadata that allows it
# to be compiled into nanopolish
import argparse
import sys
import os
from operator import itemgetter

def write_header(fh, key, value):
    fh.write("#" + key + "\t" + value + "\n")

# Argument handling
parser = argparse.ArgumentParser( description='Convert ONT model file into nanopolish format')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output-dir', type=str, required=False)
parser.add_argument('-a', '--alphabet', type=str, required=False)
args = parser.parse_args()
f = open(args.input)

# Parse metadata out of the type dir
(dirs, filename) = os.path.split(args.input)
(_, type_dir) = os.path.split(dirs)
metadata_fields = type_dir.split("_")

if(len(metadata_fields) != 5):
    sys.stderr.write("Error, could not parse type dir\n")
    sys.exit(1)

pore = metadata_fields[0]
speed = metadata_fields[2]
K = metadata_fields[3].replace("mer", "")
is_rna = type_dir.find("RNA") != -1

new_kit_name = pore + "_" + speed

alphabet = "nucleotide" if args.alphabet == "" else args.alphabet
strand = ""
if filename.find("template") != -1:
    strand = "template"
else:
    assert(filename.find("complement") != -1)
    if filename.find("pop1") != -1:
        strand = "complement.pop1"
    else:
        assert(filename.find("pop2") != -1)
        strand = "complement.pop2"


dir_str = ""
if args.output_dir is not None:
    dir_str = args.output_dir + "/"
out_name = "%s%s.%s.%smer.%s.model" % (dir_str, new_kit_name, alphabet, K, strand)

out_file = open(out_name, "w")
write_header(out_file, "ont_model_name", type_dir)
write_header(out_file, "kit", new_kit_name)
write_header(out_file, "strand", strand)
write_header(out_file, "k", K)
if args.alphabet:
    write_header(out_file, "alphabet", args.alphabet)
write_header(out_file, "original_file", type_dir + "/" + filename)

# Read k-mer states into list
states = list()

# Copy everything to the output
header = f.readline()
out_file.write(header)
for line in f:
    # ONT files shouldnt have header tags
    assert(line[0] != "#")

    fields = line.rstrip().split()

    # The ONT RNA model is in the sequencing direction, which is 5'->3'
    # Nanopolish's internal convention is to do everything 5'->3' so we reverse
    # the direction of each state here
    if is_rna:
        fields[0] = fields[0][::-1]

    states.append(fields)

for record in sorted(states, key=itemgetter(0), reverse=False):
    out_file.write("\t".join(record) + "\n")

sys.stdout.write(out_name + "\n")


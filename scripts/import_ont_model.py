#! /usr/bin/python

# This script takes a .model file provided by ONT and adds metadata that allows it
# to be compiled into nanopolish
import argparse
import sys
import os

def write_header(fh, key, value):
    fh.write("#" + key + "\t" + value + "\n")

# Argument handling
parser = argparse.ArgumentParser( description='Convert ONT model file into nanopolish format')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output-dir', type=str, required=False)
args = parser.parse_args()
f = open(args.input)

# Parse metadata out of the type dir
(dirs, filename) = os.path.split(args.input)
(_, type_dir) = os.path.split(dirs)
metadata_fields = type_dir.split("_")
if(len(metadata_fields) != 4):
    sys.stderr.write("Error, could not parse type dir\n")
    sys.exit(1)

pore = metadata_fields[0]
speed = metadata_fields[2]
K = metadata_fields[3].replace("mer", "")

new_kit_name = pore + "_" + speed

alphabet = "nucleotide"
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
write_header(out_file, "original_file", type_dir + "/" + filename)

# Copy everything to the output
for line in f:
    # ONT files shouldnt have header tags
    assert(line[0] != "#")
    out_file.write(line)

sys.stdout.write(out_name + "\n")


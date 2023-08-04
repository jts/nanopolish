#! /usr/bin/env python

import sys
import csv
import argparse

def make_key(c, s, e):
    return c + ":" + str(s) + "-" + str(e)

def parse_key(key):
    return key.replace("-", ":").split(":")

class MethylationStats:
    def __init__(self, num_reads, num_methylated, atype):
        self.num_reads = num_reads
        self.num_methylated_reads = num_methylated
        self.analysis_type = atype

    def accumulate(self, num_reads, num_methylated):
        self.num_reads += num_reads
        self.num_methylated_reads += num_methylated

    def methylation_frequency(self):
        return float(self.num_methylated_reads) / self.num_reads

def update_stats(collection, key, num_reads, num_methylated_reads, atype):
    if key not in collection:
        collection[key] = MethylationStats(num_reads, num_methylated_reads, atype)
    else:
        collection[key].accumulate(num_reads, num_methylated_reads)

def load_nanopolish(filename):
    out = dict()
    csv_reader = csv.DictReader(open(filename), delimiter='\t')

    for record in csv_reader:
        key = make_key(record["chromosome"], record["start"], record["end"])

        # skip non-singleton, for now
        if int(record["num_cpgs_in_group"]) > 1:
            continue

        num_reads = int(record["called_sites"])
        methylated_reads = int(record["called_sites_methylated"])
        out[key] = MethylationStats(num_reads, methylated_reads, "nanopolish")

    return out

def load_bisulfite(filename):
    out = dict()
    fh = open(filename)
    for line in fh:
        fields = line.rstrip().split()
        chromosome = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        num_reads = float(fields[9])
        percent_methylated = float(fields[10])
        methylated_reads = int( (percent_methylated / 100) * num_reads)
        key = ""

        # accumulate on forward strand
        if strand == "+":
            key = make_key(chromosome, str(start), str(start))
        else:
            key = make_key(chromosome, str(start - 1), str(start - 1))

        update_stats(out, key, num_reads, methylated_reads, "bisulfite")

    return out

# Load the file of methylation frequency based on the filename
def load_methylation(filename):
    if filename.find("bisulfite") != -1:
        return load_bisulfite(filename)
    else:
        return load_nanopolish(filename)

set1 = load_methylation(sys.argv[1])
set2 = load_methylation(sys.argv[2])

output = 0
print("key\tdepth_1\tfrequency_1\tdepth_2\tfrequency_2")
for key in set1:
    if key in set2:

        d1 = set1[key]
        d2 = set2[key]

        if d1.num_reads == 0 or d2.num_reads == 0:
            continue
        print("%s\t%d\t%.4f\t%d\t%.4f" % (key, d1.num_reads, d1.methylation_frequency(), d2.num_reads, d2.methylation_frequency()))
        output += 1

sys.stderr.write("set1 sites: %d set2 sites: %d output: %d\n" % (len(set1), len(set2), output))

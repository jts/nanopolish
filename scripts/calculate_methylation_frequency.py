#! /usr/bin/env python

import math
import sys
import csv
import argparse
from collections import namedtuple

class SiteStats:
    def __init__(self):
        self.num_reads = 0
        self.posterior_methylated = 0
        self.called_sites = 0
        self.called_sites_methylated = 0

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=0)
parser.add_argument('-i', '--input', type=str, required=False)
args = parser.parse_args()
assert(args.call_threshold is not None)

sites = dict()

if args.input:
    in_fh = open(args.input)
else:
    in_fh = sys.stdin
csv_reader = csv.DictReader(in_fh, delimiter='\t')

for record in csv_reader:
    num_sites = int(record['num_cpgs']) 
    # skip non-singletons for now
    if num_sites > 1:
        continue
    
    key = record['chromosome'] + ":" + record['start'] + "-" + record['end']
    if key not in sites:
        sites[key] = SiteStats()

    llr = float(record['log_lik_ratio'])

    # is the evidence strong enough at this site to make a call?
    if abs(llr) >= args.call_threshold:

        sites[key].num_reads += 1
        sites[key].called_sites += num_sites
        if llr > 0:
            sites[key].called_sites_methylated += num_sites

# header
print "\t".join(["key", "called_sites", "called_sites_methylated", "methylated_frequency"])

for key in sites:
    if sites[key].called_sites > 0:
        f = float(sites[key].called_sites_methylated) / sites[key].called_sites
        print "\t".join([str(x) for x in [key, sites[key].called_sites, sites[key].called_sites_methylated, f]])


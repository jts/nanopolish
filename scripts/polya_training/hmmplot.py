"""
Plot a random segmentation from a dataset.

Usage:
  $ python polya.out.tsv reads.fastq.readdb.index
"""
import h5py
import pandas as pd
import numpy as np
import argparse
import os
from random import choice
from collections import OrderedDict

# plotting libraries:
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns


def load_fast5_signal(read_path):
    """Load a fast5 signal from read path; return as numpy array."""
    read_h5 = h5py.File(read_path, 'r')

    # get scaling parameters:
    offset = read_h5['UniqueGlobalKey']['channel_id'].attrs['offset']
    digitisation = read_h5['UniqueGlobalKey']['channel_id'].attrs['digitisation']
    read_range = read_h5['UniqueGlobalKey']['channel_id'].attrs['range']

    # get raw integer-encoded signal:
    rn = list(read_h5['Raw']['Reads'].keys())[0]
    signal = (read_range / digitisation) * (np.array(read_h5['Raw']['Reads'][rn]['Signal']) + offset)

    # close hdf object and return numpy signal:
    read_h5.close()
    return signal


def get_state_names(header):
    """Return a list of state-start columns in the header. E.g., `[leader_start, adapter_start, ..., transcript_start]`."""
    return list(filter(lambda name: (name[-6:] == '_start'), header))


def generate_color_palette(num_colors):
    """Generate a list (of length `num_colors`) of color IDs for matplotlib."""
    # TODO(this is a hack-ish solution. Generate it mathematically!!)
    colors = ['cyan','yellow','red','green','blue', 'orange', 'green']
    return colors[:num_colors]


def main(args):
    """Filter-in PASS-ing segmentations and plot a random segmented read to file."""
    # load dataframes:
    polya = pd.read_csv(args.polya_tsv, sep='\t')
    readdb = pd.read_csv(args.readdb, sep='\t', header=None, names=['readname','location'])

    # get the names of all state-index columns:
    state_starts = get_state_names(polya.columns.values.tolist())

    # get a random read, its segmentation, and its location:
    if (args.read is None):
        row_values  = choice(polya[polya['qc_tag'] == 'PASS'][['readname', *state_starts]].values).tolist()
        read_id = row_values.pop(0)
        state_start_indices = OrderedDict()
        for k in range(len(state_starts)):
            state_start_indices[state_starts[k]] = row_values[k]
        read_path = readdb[readdb['readname'] == read_id].values[0][1]
    else:
        try:
            read_df = polya[polya['readname'] == args.read]
            row_values = choice(read_df[read_df['qc_tag'] == 'PASS'][['readname', *state_starts]].values).tolist()
            read_id = row_values.pop(0)
            state_start_indices = OrderedDict()
            for k in range(len(state_starts)):
                state_start_indices[state_starts[k]] = row_values[k]
            read_path = readdb[readdb['readname'] == read_id].values[0][1]
        except:
            raise Exception("[hmmplot.py] read id={} could not be resolved".format(args.read))

    # load fast5 file:
    signal = load_fast5_signal(read_path)

    # create dictionary of start-stop indices for each region:
    start_stop_indices = {}
    stop_idxs = [state_start_indices[name] for name in state_starts[1:]] + [signal.shape[0]]
    colors = generate_color_palette(len(state_start_indices))
    for n, (name, start_idx) in enumerate(state_start_indices.items()):
        start_stop_indices[name] = ( start_idx, stop_idxs[n], colors[n] )

    # make segmentation plot:
    plt.figure(figsize=(18,6))
    plt.plot(signal)
    for k, v in start_stop_indices.items():
        plt.axvspan(v[0], v[1], color=v[2], alpha=0.35, label=k[:-6])
    plt.legend(loc='best')
    plt.xlim(0, signal.shape[0])
    plt.title("Segmentation: {}".format(read_id))
    plt.xlabel("Sample Index (3' to 5')")
    plt.ylabel("Current (pA)")
    if (args.out is None):
        plt.savefig("segmentation.{}.png".format(read_id))
    else:
        plt.savefig(args.out)
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a random passing segmentation from a polya output file.")
    parser.add_argument("polya_tsv", help="Output TSV of `nanopolish polya {...}`")
    parser.add_argument("readdb", help="ReadDB index file from `nanopolish index {...}`")
    parser.add_argument("--out", default=None, help="Where to put the output file. [./segmentation.<READ_ID>.png]")
    parser.add_argument("--read", default=None, help="Visualize a specific read. [random read]")
    args = parser.parse_args()
    assert(os.path.exists(args.polya_tsv)), "[ERR] {} does not exist".format(args.polya_tsv)
    assert(os.path.exists(args.readdb)), "[ERR] {} does not exist".format(args.readdb)
    main(args)

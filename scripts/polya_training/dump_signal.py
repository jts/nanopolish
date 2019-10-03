"""
dump_signal.py: take a **verbose** polya-call file and dump an HDF5 file with the signal data.

A verbose polya output file can be generated with:
``$ nanopolish polya -v --reads={...} --bam={...} --genome={...} > polya.verbose.tsv``

This script outputs an HDF5 file with top-level keys given by the initial 6 characters of a read name,
under the /scaled/ group.
"""
import argparse
import os
import sys
import numpy as np
import h5py
from tqdm import tqdm, trange


# ===== Helper fns/classes =====
_DEFAULT_POLYA_HEADER = [ 'flag', 'readname', 'contig', 'pos',
                          'raw', 'scaled', 'llk_s', 'llk_l', 'llk_a',
                          'llk_p', 'llk_c', 'llk_t', 'state' ]
class PolyaIterator(object):
    """A fast, O(1)-memory class for scanning through a verbose polya output file (generated as above)."""
    def __init__(self, fname, names=_DEFAULT_POLYA_HEADER):
        """
        Args:
        * fname: path to polya file; generated with the `nanopolish polya -v` command.
        * names: [str]; the column names of the 'polya-samples' rows in the Poly(A) calls file.
        """
        self.fname = fname
        self.fhandle = open(fname, 'r')
        self.names = names

    def next(self):
        """
        Process a row of columnar entries from the poly(A) call file into a dictionary
        of values.
        """
        try:
            # get entries:
            entries = self.fhandle.readline().rstrip().split()
            # if not a polya-samples row, skip it:
            if (entries[0] != 'polya-samples'):
                return next(self)
            # format row and return:
            else:
                row = { self.names[k]: entries[k] for k in range(len(self.names)) }
                return row
        except:
            raise StopIteration()

    def close(self):
        self.fhandle.close()

    def __iter__(self):
        """Defines this as a streaming iterator."""
        return self

    def __next__(self):
        """For Python3 compatibility."""
        return self.next()

    def __del__(self):
        """Calls self.close() upon garbage collection."""
        self.close()


# ===== HDF5 interface =====
def dump_signal_hdf(args):
    """Dump an HDF5 file containing all scales samples in a verbose polyA call file."""
    # construct & open output HDF5:
    outfile = args.out if (args.out is not None) else "./samples.hdf5"
    hdf = h5py.File(outfile, 'w-') # (throw error if file already exists)
    scaled_gp = hdf.create_group('scaled')
    if args.segmentation:
        states_gp = hdf.create_group('states')

    # loop thru polya calls output file and append samples to HDF5:
    curr_read = None
    curr_samples = []
    if args.segmentation:
        curr_states = []
    for row in tqdm(PolyaIterator(args.polya)):
        # create a new read dataset based on current samples if detect a switch:
        if row['readname'] != curr_read:
            if curr_read is not None:
                try:
                    scaled_gp.create_dataset(curr_read, data=np.array(curr_samples, dtype=np.float32))
                    if args.segmentation:
                        states_gp.create_dataset(curr_read, data=np.array(curr_states, dtype='S10'))
                except:
                    pass
            # reset current read & samples
            curr_read = row['readname']
            curr_samples = []
            if args.segmentation:
                curr_states = []
            hdf.flush()
        # otherwise append raw sample:
        curr_samples.append(float(row['scaled']))
        if args.segmentation:
            curr_states.append(row['state'])
    # append final read & close HDF5 file handle:
    try:
        scaled_gp.create_dataset(curr_read, data=np.array(curr_samples, dtype=np.float32))
        if args.segmentation:
            states_gp.create_dataset(curr_read, data=np.array(curr_states, dtype='S10'))
    except:
        pass
    hdf.flush()
    hdf.close()

    # print finishing message:
    print("[dump_signal.py] HDF5 file of (scaled) picoampere signals written to: {}".format(outfile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Dump the raw signal for a set of reads.")
    parser.add_argument("polya",
                        help="A verbose `nanopolish polya` output file, containing scaled samples.")
    parser.add_argument("--out", default=None,
                        help="Output filename. [samples.hdf5]")
    parser.add_argument("--segmentation", action="store_true", default=False,
                        help="If this flag is passed, store the segmentation states as well as the scaled samples. [False]")
    args = parser.parse_args()
    # validate input args:
    assert os.path.exists(args.polya), "[dump_signal.py] File does not exist: {}".format(args.polya)
    # run main function:
    dump_signal_hdf(args)

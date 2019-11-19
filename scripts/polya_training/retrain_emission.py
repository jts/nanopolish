#! /usr/bin/env python3
"""
retrain_emission.py: take an HDF5 file and segmentations, and output parameters of a mixture model.
"""
# std lib:
import argparse
import os
import sys
import random
from collections import defaultdict
from tqdm import tqdm
# numerics:
import numpy as np
import h5py
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture


def pool_reads(h, K):
    """
    Select (up to) K random segmented reads from the dataset `h`. Return as a dictionary of pooled scaled samples
    of form { <REGION_NAME> :: str -> <SAMPLES> :: NDArray(float) }.
    """
    # collect scaled samples for each state:
    pool = defaultdict(list)
    rnames = random.sample(h['scaled'].keys(), min(K, len(h['scaled'].keys())))
    for rid in tqdm(rnames):
        try:
            assert(len(h['scaled'][rid]) == len(h['states'][rid]))
            for k in range(len(h['states'][rid])):
                pool[ h['states'][rid][k] ].append( h['scaled'][rid][k] )
        except:
            pass

    # process into a dict of numpy arrays and return:
    pool = dict(pool)
    for k, v in pool.items():
        pool[k] = np.array(v)
    return pool


def retrain_emission(hdf_path, nreads, bayesian, components, verbose):
    """Retrain gaussian mixture model from parameters."""
    # load dataset:
    hdf = h5py.File(args.hdf_path, 'r')
    assert ('states' in hdf.keys() and 'scaled' in hdf.keys()), \
        "[retrain_emission.py] ERR: both `samples` and `states` must be groups in the HDF5."

    # select up to `nreads` random segmented reads from the dataset:
    print("[retrain_emission.py] Collecting and pooling {} random reads (this may take a while...)".format(nreads))
    segments = pool_reads(hdf, nreads)

    # compute GMM parameters for each segment:
    CONFIG = {
        'ncomp': components,
        'niter': 100,
        'ninit': 5,
        'verbose': (1 if verbose else 0),
        'bayesian': bayesian
    }
    print("----- TRAINING CONFIG -----")
    for k,v in CONFIG.items():
        print("* {0} = {1}".format(k,v))
    gmm = {}
    for k,v in segments.items():
        if v.shape[0] < 10:
            print("[retrain_emissions.py] Fewer than 10 samples for state {}; skipping...".format(k))
            pass
        # train GMM:
        if CONFIG['bayesian']:
            gmm[k] = BayesianGaussianMixture(
                n_components=CONFIG['ncomp'], max_iter=CONFIG['niter'], n_init=CONFIG['ninit'],
                verbose=CONFIG['verbose']).fit(v.reshape(-1,1))
        else:
            gmm[k] = GaussianMixture(
                n_components=CONFIG['ncomp'], max_iter=CONFIG['niter'], n_init=CONFIG['ninit'],
                verbose=CONFIG['verbose']).fit(v.reshape(-1,1))

    # print mixture model properties for each segment:
    for k,v in gmm.items():
        print("===== [{}] =====".format(k))
        print("* Weights: {}".format(v.weights_))
        print("* Means: {}".format(v.means_))
        print("* Covariances: {}".format(v.covariances_))

    hdf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Train a mixture model.")
    parser.add_argument("hdf_path",
                        help="Path to HDF5 file with segmented signal paths.")
    parser.add_argument("--nreads", default=50, type=int,
                        help="Number of random reads to pool together and retrain upon. [50]")
    parser.add_argument("--bayesian", default=False, action='store_true',
                        help="Use a dirichlet process mixture model. [False]")
    parser.add_argument("--verbose", default=False, action='store_true',
                        help="Print verbose outputs during training. [False]")
    parser.add_argument("--components", default=2, type=int,
                        help="If DPMM, max components; else fixed number of GMM components. [2]")
    args = parser.parse_args()
    assert (os.path.exists(args.hdf_path)), "File does not exist: {}".format(args.hdf_path)
    retrain_emission(args.hdf_path, args.nreads, args.bayesian, args.components, args.verbose)

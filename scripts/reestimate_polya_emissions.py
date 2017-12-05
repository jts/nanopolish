"""
reestimate_polya_emissions.py: given two `polya-samples` TSV files based on different
underlying kmer models (with the newer TSV giving failing poly(A) segmentations),
infer the best new parameters for the HMM emissions.

Usage:
$ python reestimate_polya_emissions.py samples.old.tsv seg.old.tsv samples.new.tsv
where:
* `samples.old.tsv` is the output of `nanopolish polya -vv [...] | grep 'polya-samples'`,
with the **old** kmer models;
* `seg.old.tsv` is the output of `nanopolish polya -v [...] | grep 'polya-segmentation'`,
with the **old** kmer models;
* `samples.new.tsv` is the output of `nanopolish polya -vv [...] | grep 'polya-samples'`,
with the **new** kmer models.

Dependencies:
* numpy >= 1.11.2
* scipy >= 0.18.1
"""
import csv
import numpy as np
import argparse
import os
from scipy.stats import norm


log_inv_sqrt_2pi = np.log(0.3989422804014327)
def log_normal_pdf(xs, mu, sigma):
    """Compute the log-normal PDF of a given sample(s) against a mu and sigma."""
    alpha = (xs - mu) * np.reciprocal(sigma)
    return ( log_inv_sqrt_2pi - np.log(sigma) + (-0.5 * alpha * alpha) )


def fit_gaussian(samples):
    """Given a numpy array of floating point samples, fit a gaussian distribution."""
    mu, sigma = norm.fit(samples)
    return (mu,sigma)


def old_tsv_to_numpy(tsv_path):
    """
    Read a TSV containing raw samples and return a dictionary consisting
    of the following numpy datasets:
    * L_loglkhd: the log-likelihoods of the samples belonging to the LEADER segment.
    * A_loglkhd: the log-likelihoods of the samples belonging to the ADAPTER segment.
    * P_loglkhd: the log-likelihoods of the samples belonging to the POLYA segment.
    * T_loglkhd: the log-likelihoods of the samples belonging to the TRANSCRIPT segment.
    """
    # instantiate arrays to hold values:
    L_loglkhd = []
    A_loglkhd = []
    P_loglkhd = []
    T_loglkhd = []

    # loop over TSV file and append data to arrays:
    str2int = { 'LEADER': 0, 'ADAPTER': 1, 'POLYA': 2, 'TRANSCRIPT': 3 }
    with open(tsv_path, 'r') as f:
        headers =  ['tag','read_id', 'chr', 'idx', 'sample', 'scaled_sample',
                    'l_llh','a_llh','p_llh','t_llh','region']
        rdr = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE, fieldnames=headers)
        for row in rdr:
            # parse row fields:
            l_llh = float(row['l_llh'])
            a_llh = float(row['a_llh'])
            p_llh = float(row['p_llh'])
            t_llh = float(row['t_llh'])
            region = row['region']
            
            # append log-likelihoods to appropriate arrays:
            if region == 'LEADER':
                L_loglkhd.append(l_llh)
            if region == 'ADAPTER':
                A_loglkhd.append(a_llh)
            if region == 'POLYA':
                P_loglkhd.append(p_llh)
            if region == 'TRANSCRIPT':
                T_loglkhd.append(t_llh)
    
    return { "L_loglkhd": np.array(L_loglkhd, dtype=float),
             "A_loglkhd": np.array(A_loglkhd, dtype=float),
             "P_loglkhd": np.array(P_loglkhd, dtype=float),
             "T_loglkhd": np.array(T_loglkhd, dtype=float) }


def make_segmentation_dict(segmentations_tsv_path):
    """
    Load a segmentations TSV file. Rows of `segmentations_tsv_path` look like this:

    tag                 read_id: pos:       A_0:    P_0:     P_1:     RR:    P(A)L:
    polya-segmentation  fc06...  161684804  1851.0  8354.0   11424.0  73.76  75.18

    Note that this function only takes the first available segmentation for each read, i.e.
    if a read id appears more than once in the TSV, only the first segmentation is kept, and
    later occurrences of the read id in the TSV are ignored.
    """
    segments = {}
    # loop thru TSV and update the list of segmentations:
    with open(segmentations_tsv_path, 'r') as f:
        headers = ['tag', 'read_id', 'pos', 'A_start', 'P_start', 'P_end', 'rate', 'length']
        rdr = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE, fieldnames=headers)
        for row in rdr:
            if row['read_id'] not in segments.keys():
                segments[row['read_id']] = { 'A_start': int(float(row['A_start'])),
                                             'P_start': int(float(row['P_start'])),
                                             'P_end': int(float(row['P_end'])) }
    return segments


def region_search(read_id, sample_ix, segmentations):
    """
    Given a dictionary of ("gold-standard") segmentations, look up the region that a
    given read and sample index belongs to.

    Returns an integer label out of 0,1,2,3,4, where:
    0 => LEADER, 1 => ADAPTER, 2 => POLYA, 3 => TRANSCRIPT, 4 => UNKNOWN
    """
    # find read ID in segmentations:
    read_key = None
    for long_read_id in segmentations.keys():
        if long_read_id[0:len(read_id)] == read_id:
            read_key = long_read_id

    # return UNK if read not found:
    if read_key == None:
        return 4

    # find region that `sample_ix` belongs to:
    a_start = segmentations[read_key]['A_start']
    p_start = segmentations[read_key]['P_start']
    p_end = segmentations[read_key]['P_end']
    if (sample_ix < a_start):
        return 0
    if (sample_ix < p_start):
        return 1
    if (sample_ix <= p_end):
        return 2
    if (sample_ix > p_end):
        return 3
    return 4


def new_tsv_to_numpy(tsv_path, segmentations):
    """
    Read a TSV of new, miscalled samples and a dictionary of correct segmentations (coming from
    an older, correct TSV) and return a dict of numpy arrays.

    Args:
    * tsv_path: ...
    * segmentations: a dictionary of segmentation intervals, given by ...

    Returns: a dictionary of numpy arrays.
    [TODO: describe the dict/finish this multiline comment]
    """
    # instantiate arrays to hold values:
    L_samples = []
    A_samples = []
    P_samples = []
    T_samples = []

    # loop over TSV file and append data to arrays:
    with open(tsv_path, 'r') as f:
        headers =  ['tag','read_id', 'chr', 'idx', 'sample', 'scaled_sample',
                    'l_llh','a_llh','p_llh','t_llh','region']
        rdr = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE, fieldnames=headers)
        for row in rdr:
            scaled_sample = float(row['scaled_sample'])
            read = row['read_id']
            contig = row['chr']
            index = int(row['idx'])
            region = region_search(read, index, segmentations)
            if region == 0:
                L_samples.append(scaled_sample)
            if region == 1:
                A_samples.append(scaled_sample)
            if region == 2:
                P_samples.append(scaled_sample)
            if region == 3:
                T_samples.append(scaled_sample)
            
    return { "L_samples": np.array(L_samples, dtype=float),
             "A_samples": np.array(A_samples, dtype=float),
             "P_samples": np.array(P_samples, dtype=float),
             "T_samples": np.array(T_samples, dtype=float) }


def compare_likelihoods(old_log_lkhds, new_samples, new_mu, new_sigma):
    """
    Given log-likelihood data from old and new HMM settings, benchmark the accuracy of
    new log likelihoods against old one.

    Args:
    * old_log_lkhds: 1D np.array containing the old log-likelihoods for this region.
    * new_samples: 1D np.array containing the new scaled samples for this region.
    * new_mu: float representing the newly-inferred mean of this region.
    * new_sigma: float representing the newly-inferred standard deviation of this region.

    Returns: a tuple (avg_old_llh, avg_new_llh) where:
    * avg_old_llh: float representing the average of `old_log_lkhds`
    * avg_new_llh: float representing the average of the newly computed log-likelihoods based
    on a gaussian with parameters `(new_mu, new_sigma)`.
    """
    # compute log-likelihoods of new values:
    new_log_lkhds = log_normal_pdf(new_samples, new_mu, new_sigma)

    # return averaged values:
    avg_old_llh = np.mean(old_log_lkhds)
    avg_new_llh = np.mean(new_log_lkhds)
    return (avg_old_llh, avg_new_llh)


def main(old_samples_tsv, old_segmentations_tsv, new_samples_tsv, benchmark=True):
    """
    Infer and print the new values for mu and sigma (for each of L, A, P, T) to STDOUT.

    Args:
    * old_samples_tsv: path to TSV file containing polya-samples data from an older kmer model.
    * old_segmentations_tsv: path to TSV file containing polya-segmentation data from an older kmer model.
    * new_samples_tsv: path to TSV file containing polya-samples data from the newer kmer model.

    Returns: N/A, prints outputs to STDOUT.
    """
    ### read all samples into numpy arrays:
    print("Loading data from TSV...")
    old_data = old_tsv_to_numpy(old_samples_tsv)
    segmentations = make_segmentation_dict(old_segmentations_tsv)
    new_data = new_tsv_to_numpy(new_samples_tsv, segmentations)
    print("... Datasets loaded.")

    ### infer best possible new mu,sigma for each of L, A, P, T:
    print("Fitting gaussians to new scaled samples (this may take a while)...")
    new_mu_L, new_sigma_L = fit_gaussian(new_data['L_samples'])
    new_mu_A, new_sigma_A = fit_gaussian(new_data['A_samples'])
    new_mu_P, new_sigma_P = fit_gaussian(new_data['P_samples'])
    new_mu_T, new_sigma_T = fit_gaussian(new_data['T_samples'])

    ### print to stdout:
    print("New params for LEADER: mu = {0}, var = {1}, stdv = {2}".format(new_mu_L, new_sigma_L, np.sqrt(new_sigma_L)))
    print("New params for ADAPTER: mu = {0}, var = {1}, stdv = {2}".format(new_mu_A, new_sigma_A, np.sqrt(new_sigma_A)))
    print("New params for POLYA: mu = {0}, var = {1}, stdv = {2}".format(new_mu_P, new_sigma_P, np.sqrt(new_sigma_P)))
    print("New params for TRANSCRIPT: mu = {0}, var = {1}, stdv = {2}".format(new_mu_T, new_sigma_T, np.sqrt(new_sigma_T)))

    ### optionally, benchmark:
    if not benchmark:
        return

    print("===== Emission Log-Likelihood Benchmarks =====")
    old_L_llh, new_L_llh = compare_likelihoods(old_data['L_loglkhd'], new_data['L_samples'], new_mu_L, new_sigma_L)
    print("> Average LEADER log-probs:")
    print("> Old avg. log-likelihood: {0} | New avg. log-likelihood: {1}".format(old_L_llh, new_L_llh))

    old_A_llh, new_A_llh = compare_likelihoods(old_data['A_loglkhd'], new_data['A_samples'], new_mu_A, new_sigma_A)
    print("> Average ADAPTER log-probs:")
    print("> Old avg. log-likelihood: {0} | New avg. log-likelihood: {1}".format(old_A_llh, new_A_llh))

    old_P_llh, new_P_llh = compare_likelihoods(old_data['P_loglkhd'], new_data['P_samples'], new_mu_P, new_sigma_P)
    print("> Average POLYA log-probs:")
    print("> Old avg. log-likelihood: {0} | New avg. log-likelihood: {1}".format(old_P_llh, new_P_llh))

    old_T_llh, new_T_llh = compare_likelihoods(old_data['T_loglkhd'], new_data['T_samples'], new_mu_T, new_sigma_T)
    print("> Average TRANSCRIPT log-probs:")
    print("> Old avg. log-likelihood: {0} | New avg. log-likelihood: {1}".format(old_T_llh, new_T_llh))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Infer new Poly(A) emission parameters.")
    parser.add_argument("old_samples_tsv", help="Path to TSV file of old samples.")
    parser.add_argument("segmentation_tsv", help="Path to segmentations for reads.")
    parser.add_argument("new_samples_tsv", help="Path to TSV file of new samples.")
    parser.add_argument("--benchmark", default=True, type=bool, dest="benchmark",
                        help="If `--benchmark=False`, don't the new estimated HMM parameters.")
    args = parser.parse_args()
    # sanity checks:
    assert os.path.exists(args.old_samples_tsv)
    assert os.path.exists(args.segmentation_tsv)
    assert os.path.exists(args.new_samples_tsv)
    # run inference and (optional) benchmarking of new parameters:
    main(args.old_samples_tsv, args.segmentation_tsv, args.new_samples_tsv, benchmark=args.benchmark)

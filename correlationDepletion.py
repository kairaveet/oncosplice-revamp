import pandas as pd
import numpy as np


def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)


def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)


def weighted_corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))


def calculate_sample_weights(full_psi_file):

    percentage_non_missing = (1 - full_psi_file.isna().mean()) * 100

    # The mean() method computes the mean of the boolean mask for each column, which gives the proportion of NaN values.

    weights = pd.DataFrame(percentage_non_missing, columns=['Percentage non-missing values'])

    return weights


def calculate_weighted_correlation(nmf_basis_matrix, full_psi_file, weights):

    nmf_basis_matrix = pd.DataFrame(nmf_basis_matrix).transpose()  # the shape of this matrix should be number of ranks by number of samples
    nmf_basis_matrix = nmf_basis_matrix.set_axis(full_psi_file.columns, axis=1)  # make sure the column names of the nmf matrix is the same as sample names (columns of psi file)

    weights = np.array(weights.loc[:, 'Percentage non-missing values'])

    corr_matrix = np.zeros([np.shape(nmf_basis_matrix)[0], np.shape(full_psi_file)[0]])  # create an empty array (of zeros) of size number of clusters by number of splicing events

    for cluster in range(np.shape(nmf_basis_matrix)[0]):
        for event in range(np.shape(full_psi_file)[0]):
            corr_matrix[cluster, event] = weighted_corr(nmf_basis_matrix.iloc[cluster, :], full_psi_file.iloc[event, :], weights)  # this is a slow step -- how can you optimize this?

    return corr_matrix


def deplete_events(nmf_basis_matrix, full_psi_file, corr_threshold=0.3, strictness="easy"):

    weights = calculate_sample_weights(full_psi_file)

    corr_matrix = calculate_weighted_correlation(nmf_basis_matrix, full_psi_file, weights)

    corr_matrix = np.array(corr_matrix)

    mask = abs(corr_matrix) < corr_threshold

    if strictness == "easy":
        has_low_corr = np.all(mask, axis=0)
    else:
        has_low_corr = np.any(mask, axis=0)  # should this check for the TRUE in at least one cluster or in ALL clusters?

    depleted_psi_file = full_psi_file.loc[~has_low_corr]

    return depleted_psi_file

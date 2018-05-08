# coding=utf-8

import pandas as pd
from sklearn import mixture
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt


__alls__ = ['binarize']


def _derive_threshold(auc_mtx, regulon_name):
    # We fit a Gaussian mixture of 2 components on the AUC distribution using EM algorithm. The threshold is defined
    # as the largest mean minus the standard deviation of the guassian that corresponds to this Gaussian component.
    data = auc_mtx[regulon_name].values.reshape(-1, 1)
    gmm = mixture.GaussianMixture(n_components=2, covariance_type='full').fit(data)
    idx = np.argmax(gmm.means_)
    return max((gmm.means_[idx] - 1 * sqrt(gmm.covariances_[idx]))[0], 0)

def binarize(auc_mtx: pd.DataFrame) -> (pd.DataFrame, pd.Series):
    """

    "Binarize" the supplied AUC matrix, i.e. decide if for each cells in the matrix a regulon is active or not based
    on the bimodal distribution of the AUC values for that regulon.

    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :return: A "binarized" dataframe and a series containing the AUC threshold used for each regulon.
    """
    def derive_thresholds(auc_mtx):
        return pd.Series(index=auc_mtx.columns, data=[_derive_threshold(auc_mtx, name) for name in auc_mtx.columns])
    thresholds = derive_thresholds(auc_mtx)
    return (auc_mtx > thresholds).astype(int), thresholds

def plot_binarization(auc_mtx, regulon_name, bins=200):
    """

    Plot the "binarization" process for the given regulon.

    :param auc_mtx:
    :param regulon_name:
    :param bins:
    :return:
    """
    auc_mtx[regulon_name].hist(bins=bins)
    threshold = _derive_threshold(auc_mtx, regulon_name)
    ylim = plt.ylim()
    plt.plot([threshold]*2, ylim, 'r:')
    plt.ylim(ylim)
    plt.xlabel('AUC')
    plt.ylabel('#')
    plt.title(regulon_name)

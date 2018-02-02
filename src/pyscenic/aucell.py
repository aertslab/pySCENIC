# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from itertools import repeat
from .genesig import Regulome


def create_rankings(ex_mtx: pd.DataFrame) -> pd.DataFrame:
    """
    Create a whole genome rankings dataframe from a single cell expression profile dataframe.

    :param ex_mtx: The expression profile matrix. The columns should correspond to different cells, the rows to different
        genes.
    :return: A genome rankings dataframe.
    """
    return ex_mtx.rank(axis=0, ascending=False, method='first').astype('int64')


def enrichment(rnk_mtx: pd.DataFrame, regulome: Regulome, rank_threshold: int = 5000, auc_threshold: float = 0.05) -> pd.DataFrame:
    """
    Calculate the enrichment of the regulome for the cells in the ranking dataframe.

    :param rnk_mtx:
    :param regulome:
    :param rank_threshold:
    :param auc_threshold:
    :return:
    """

    # Load rank of genes from database.
    total_genes = len(rnk_mtx.columns)
    rank_cutoff = int(round(auc_threshold * total_genes))
    df = rnk_mtx[rnk_mtx.index.isin(regulome.genes)]
    features, genes, rankings = df.columns.values, df.index.values, df.T.values
    weights = np.asarray([regulome[gene] for gene in genes] + [0.0])

    # Calculate recovery curves.
    def calc_rcc(ranking, weights, total_genes, rank_threshold):
        curranking = np.append(ranking, total_genes)
        return np.cumsum(np.bincount(curranking, weights=weights)[:rank_threshold])
    # Apply along axis does not improve performance, only more readable code.
    rccs = np.apply_along_axis(calc_rcc, 1, rankings, weights, total_genes, rank_threshold)

    # Calculate AUC.
    maxauc = float(rank_cutoff * total_genes)
    aucs = rccs[:, :rank_cutoff].sum(axis=1) / maxauc

    return pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(features, repeat(regulome.transcription_factor))),
                                                        names=["Cell", "Regulome"]),
                        data={"AUC": aucs})
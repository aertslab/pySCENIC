# -*- coding: utf-8 -*-

import numpy
import pandas as pd
from itertools import repeat

from .rnkdb import RankingDatabase
from .genesig import GeneSignature


def enrichment(rnkdb: RankingDatabase, gs: GeneSignature, rank_threshold: int = 5000, auc_threshold: float = 0.05) -> pd.DataFrame:
    """
    Calculate AUC and NES for all regulatory features in the supplied database using the genes of the give signature.

    :param rnkdb: The database.
    :param gs: The gene signature.
    :param rank_threshold:
    :param auc_threshold:
    :return: A dataframe containing all information.
    """
    assert rnkdb, "A database must be supplied"
    assert gs, "A gene signature must be supplied"
    assert 0 < rank_threshold < rnkdb.total_genes, \
            "Rank threshold must be an integer between 1 and {0:d}".format(rnkdb.total_genes)
    assert 0.0 < auc_threshold <= 1.0, "AUC threshold must be a fraction between 0.0 and 1.0"
    rank_cutoff = int(round(auc_threshold * rnkdb.total_genes))
    assert rank_cutoff <= rank_threshold, \
            "An AUC threshold of {0:f} corresponds to {1:d} top ranked genes/regions in the database. " \
            "Please increase the rank threshold or decrease the AUC threshold.".format(auc_threshold, rank_cutoff)

    # Load rank of genes from database.
    features, genes, weights, rankings = rnkdb.load(gs)

    # Calculate recovery curves.
    def calc_rcc(ranking, weights, total_genes, rank_threshold):
        curranking = numpy.append(ranking, total_genes)
        return numpy.cumsum(numpy.bincount(curranking, weights=weights)[:rank_threshold])
    # Apply along axis does not improve performance, only more readable code.
    rccs = numpy.apply_along_axis(calc_rcc, 1, rankings, weights, rnkdb.total_genes, rank_threshold)

    # Calculate AUC and NES.
    maxauc = float(rank_cutoff * rnkdb.total_genes)
    aucs = rccs[:, :rank_cutoff].sum(axis=1) / maxauc
    ness = (aucs - aucs.mean()) / aucs.std()

    # Calculate average recovery curve.
    #avgrcc = numpy.average(rccs, axis=1)
    #stdrcc = numpy.std(rccs, axis=1)
    #avg2stdrcc = avgrcc + 2.0 * stdrcc

    df_nes = pd.DataFrame(index=features, columns=[("Enrichment", "AUC"), ("Enrichment", "NES")], data=[aucs, ness])
    df_rnks = pd.DataFrame(index=features, columns=list(zip(repeat("Ranking"), genes)), data=rankings)
    df_rccs = pd.DataFrame(index=features, columns=list(zip(repeat("Recovery"), numpy.arange(rank_threshold))), data=rccs)
    return pd.concat([df_nes, df_rccs, df_rnks], axis=1)


def leading_edge(row, avg2stdrcc):

    ranking = row["Ranking"]
    genes =
    rcc = row[""]
    #TODO:
    # Use partial function application to make this function really appliable to the rows of a dataframe.
    def critical_point(rcc, avg2stdrcc, rank_threshold):
        """ Returns (x,y). """
        x_values = numpy.arange(1, rank_threshold + 1)
        y_values = rcc - avg2stdrcc
        y_max = y_values.max()
        x_max = int(x_values[y_values == y_max][0])
        return x_max, rcc[x_max - 1]

    def get_genes(genes, ranking, rank):
        sorted_idx = numpy.argsort(ranking)
        ranking = ranking[sorted_idx]
        gene_ids = genes[sorted_idx]
        filtered_idx = ranking < rank
        return zip(ranking[filtered_idx] + 1, gene_ids[filtered_idx])

    rank_threshold = ranking.shape
    rank, n_recovered_genes = critical_point(rcc, avg2stdrcc, rank_threshold)
    return get_genes(genes, ranking, rank)

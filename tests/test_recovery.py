# -*- coding: utf-8 -*-

from pyscenic.recovery import enrichment4features as enrichment, auc1d, weighted_auc1d, rcc2d

import os
import numpy as np
from configparser import ConfigParser
from pyscenic.rnkdb import SQLiteRankingDatabase as RankingDatabase
from pyscenic.genesig import GeneSignature


TEST_DATABASE = "hg19-500bp-upstream-10species"
TEST_SIGNATURE = "msigdb_cancer_c6"


def load_db_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_sqlitedb.ini'))
    return config[section]

def load_gs_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_genesig.ini'))
    return config[section]

def test_enrichment():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[0]
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    df = enrichment(db, gs)

def test_auc1d_1():
    # Check if AUC is calculated correctly when a gene is recovered at the rank threshold.
    # In the python implementation it should be included in the AUC calculation.
    # For the current gene list the non-normalized AUC should be (2*1)+(2*2)+(2*3)+(1*4) = 16.
    total_genes = 100
    auc_rank_threshold = 8 # This rank threshold is inclusive.
    ranking = np.asarray([2, 4, 6, 8]) - 1 # The databases have a zero-based ranking system
    auc_max = auc_rank_threshold * total_genes
    assert 16.0/800 == auc1d(ranking, auc_rank_threshold, auc_max)

def test_auc1d_2():
    # For the current gene list the non-normalized AUC should be (2*1)+(2*2)+(3*3) = 15.
    total_genes = 100
    auc_rank_threshold = 8 # This rank threshold is inclusive.
    ranking = np.asarray([2, 4, 6]) - 1 # The databases have a zero-based ranking system
    auc_max = auc_rank_threshold * total_genes
    assert 15.0/800 == auc1d(ranking, auc_rank_threshold, auc_max)

def test_weighted_auc1d():
    # CAVE: In python the ranking databases are 0-based. The only difference a 1-based system has on the calc
    # of the AUC is that in the latter the rank threshold would not be included. This has an infuence on the
    # normalization factor max AUC but nothing else.
    total_genes = 100
    auc_rank_threshold = 8
    ranking = np.asarray([2, 4, 6]) - 1 # The databases have a zero-based ranking system
    weights = np.ones(len(ranking))
    auc_max = 1.0 # Disable normalization.
    assert 15.0 == weighted_auc1d(ranking, weights, auc_rank_threshold, auc_max)


def test_weighted_auc1d_batch():
    # The assumption taken here is that for weights uniformely being set to 1.0, auc1d and
    # weighted_auc1d should always have the same output.
    N = 100000
    total_genes = 100
    for _ in range(N):
        auc_rank_threshold = np.random.randint(2, high=total_genes)
        n_genes = np.random.randint(20, high=total_genes)
        ranking = np.random.randint(low=0, high=total_genes, size=n_genes)
        weights = np.ones(n_genes)
        auc_max = 1.0 # Disable normalization.
        assert auc1d(ranking, auc_rank_threshold, auc_max) == weighted_auc1d(ranking, weights, auc_rank_threshold, auc_max)

def test_weighted_rcc2d_batch():
    # The assumption taken here is that for weights uniformely being set to 1.0, auc1d and
    # weighted_auc1d should always have the same output.
    N = 100000
    total_genes = 100
    for _ in range(N):
        #rccs[row_idx, :] = np.cumsum(np.bincount(curranking, weights=weights)[:rank_threshold])
        #ValueError: could not broadcast input array from shape (89) into shape (97)
        auc_rank_threshold = np.random.randint(2, high=total_genes)
        n_genes = np.random.randint(20, high=total_genes)
        ranking = np.random.randint(low=0, high=total_genes, size=n_genes)
        rankings = ranking.reshape((1, n_genes))
        rankings = np.append(rankings, np.full(shape=(1, 1), fill_value=total_genes), axis=1)
        weights = np.ones(n_genes)
        auc_max = 1.0 # Disable normalization.
        assert rcc2d(rankings, np.insert(weights, len(weights), 0.0), total_genes)[:, :auc_rank_threshold].sum(axis=1) == weighted_auc1d(ranking, weights, auc_rank_threshold, auc_max)

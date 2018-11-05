# -*- coding: utf-8 -*-

from pyscenic.recovery import enrichment4features as enrichment, auc1d, weighted_auc1d, rcc2d

import pytest
import numpy as np
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.genesig import GeneSignature
from pkg_resources import resource_filename


TEST_DATABASE_FNAME = resource_filename('resources.tests', "hg19-tss-centered-10kb-10species.mc9nr.feather")
TEST_DATABASE_NAME = "hg19-tss-centered-10kb-10species"
TEST_SIGNATURE_FNAME = resource_filename('resources.tests', "c6.all.v6.1.symbols.gmt")


@pytest.fixture
def db():
    return RankingDatabase(TEST_DATABASE_FNAME, TEST_DATABASE_NAME,)


@pytest.fixture
def gs():
    return GeneSignature.from_gmt(TEST_SIGNATURE_FNAME,
                                  gene_separator="\t", field_separator="\t", )[0]


def test_enrichment(db, gs):
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

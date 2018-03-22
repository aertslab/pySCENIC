# coding=utf-8

import pytest

import pandas as pd

from pyscenic.genesig import GeneSignature
from pyscenic.aucell import derive_auc_threshold, aucell, create_rankings
from pkg_resources import resource_filename



NOMENCLATURE = "HGNC"
TEST_EXPRESSION_MTX_FNAME = resource_filename('resources', "GSE103322.em.hgnc.sample.cxg.csv")
TEST_SIGNATURE_FNAME = resource_filename('resources', "c6.all.v6.1.symbols.gmt")


@pytest.fixture
def exp_matrix():
    return pd.read_csv(TEST_EXPRESSION_MTX_FNAME, sep=',', header=0, index_col=0)


@pytest.fixture
def gs():
    return GeneSignature.from_gmt(TEST_SIGNATURE_FNAME, NOMENCLATURE,
                                  gene_separator="\t", field_separator="\t", )


def test_create_rankings():
    ex_mtx = exp_matrix()
    df_rnk = create_rankings(ex_mtx)
    n_genes = ex_mtx.shape[1]
    assert len(df_rnk.sum(axis=1).unique()) == 1
    assert (df_rnk + 1).sum(axis=1).unique()[0] == (n_genes * (n_genes+1))/2.0


def test_aucell():
    ex_mtx = exp_matrix()
    percentiles = derive_auc_threshold(ex_mtx)
    aucs_mtx = aucell(ex_mtx, gs(), auc_threshold=percentiles[0.01], num_cores=1)




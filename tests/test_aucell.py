# coding=utf-8

import pandas as pd
import pytest
from ctxcore.genesig import GeneSignature
from pkg_resources import resource_filename

from pyscenic.aucell import aucell, create_rankings, derive_auc_threshold

NOMENCLATURE = "HGNC"
TEST_EXPRESSION_MTX_FNAME = resource_filename(
    "resources.tests", "GSE103322.em.hgnc.sample.cxg.csv"
)
TEST_SIGNATURE_FNAME = resource_filename("resources.tests", "c6.all.v6.1.symbols.gmt")


@pytest.fixture
def exp_matrix():
    return pd.read_csv(TEST_EXPRESSION_MTX_FNAME, sep=",", header=0, index_col=0)


@pytest.fixture
def gs():
    return GeneSignature.from_gmt(
        TEST_SIGNATURE_FNAME,
        gene_separator="\t",
        field_separator="\t",
    )


def test_create_rankings(exp_matrix):
    df_rnk = create_rankings(exp_matrix)
    n_genes = exp_matrix.shape[1]
    assert len(df_rnk.sum(axis=1).unique()) == 1
    assert (df_rnk + 1).sum(axis=1).unique()[0] == (n_genes * (n_genes + 1)) / 2.0


def test_aucell_w1(exp_matrix, gs):
    percentiles = derive_auc_threshold(exp_matrix)
    aucs_mtx = aucell(exp_matrix, gs, auc_threshold=percentiles[0.01], num_workers=1)


def test_aucell_w2(exp_matrix, gs):
    percentiles = derive_auc_threshold(exp_matrix)
    aucs_mtx = aucell(exp_matrix, gs, auc_threshold=percentiles[0.01], num_workers=4)


def test_aucell_mismatch(exp_matrix, gs):
    percentiles = derive_auc_threshold(exp_matrix)
    gss = [
        GeneSignature(name="test", gene2weight=list(map("FAKE{}".format, range(100))))
    ] + gs
    aucs_mtx = aucell(exp_matrix, gss, auc_threshold=percentiles[0.01], num_workers=1)
    print(aucs_mtx.head())

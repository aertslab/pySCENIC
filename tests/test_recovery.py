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

def test_auc1d():
    total_genes = 100
    auc_rank_threshold = 8
    ranking = np.asarray([2, 4, 6])
    weights = np.ones(len(ranking))
    auc_max = 1.0 # Disable normalization.
    assert 12.0 == auc1d(ranking, auc_rank_threshold, auc_max)

def test_weighted_auc1d():
    total_genes = 100
    auc_rank_threshold = 8
    ranking = np.asarray([2, 4, 6])
    weights = np.ones(len(ranking))
    auc_max = 1.0 # Disable normalization.
    assert 12.0 == weighted_auc1d(ranking, weights, auc_rank_threshold, auc_max)


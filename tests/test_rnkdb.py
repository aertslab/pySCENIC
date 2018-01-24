# -*- coding: utf-8 -*-

import pytest
import os
from configparser import ConfigParser
from pyscenic.rnkdb import RankingDatabase
from pyscenic.genesig import GeneSignature, Regulome


TEST_DATABASE = "hg19-500bp-upstream-10species"
TEST_SIGNATURE = "msigdb_cancer_c6"


def load_db_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_rnkdb.ini'))
    return config[section]

def load_gs_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_genesig.ini'))
    return config[section]

def test_init():
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    assert db.name == "hg19-500bp-upstream-10species"
    assert db.nomenclature == "HGNC"

def test_load_full():
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    features, genes, rankings = db.load_full()
    assert len(features) == 24453
    assert len(genes) == 22284
    assert rankings.shape[0] == len(features)
    assert rankings.shape[1] == len(genes)

def test_load():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[0]
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    features, genes, rankings = db.load(gs)
    assert len(features) == 24453
    assert len(genes) == 29
    assert rankings.shape[0] == len(features)
    assert rankings.shape[1] == len(genes)

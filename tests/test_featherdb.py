# -*- coding: utf-8 -*-

import pytest
import os
from configparser import ConfigParser
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase, convert2feather
from pyscenic.genesig import GeneSignature


TEST_DATABASE = "hg19-500bp-upstream-10species"
TEST_SIGNATURE = "msigdb_cancer_c6"


def load_db_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_featherdb.ini'))
    return config[section]

@pytest.fixture
def db():
    cfg = load_db_info(TEST_DATABASE)
    if not os.path.exists(cfg['fname']):
        fname = cfg['fname']
        convert2feather(fname=os.path.basename(fname),
                        name=cfg['name'],
                        out_folder="{}.db".format(os.path.dirname(fname)),
                        nomenclature=cfg['nomenclature'])
    return RankingDatabase(**cfg)

def load_gs_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_genesig.ini'))
    return config[section]

@pytest.fixture
def gs():
    return GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[0]

def test_init(db):
    assert db.name == "hg19-500bp-upstream-10species"
    assert db.nomenclature == "HGNC"

def test_total_genes(db):
    assert db.total_genes == 22284

def test_load_full(db):
    rankings = db.load_full()
    assert len(rankings.index) == 24453
    assert len(rankings.columns) == 22284

def test_load(db, gs):
    rankings = db.load(gs)
    assert len(rankings.index) == 24453
    assert len(rankings.columns) == 29

# -*- coding: utf-8 -*-

from pyscenic.recovery import enrichment

import os
from configparser import ConfigParser
from pyscenic.rnkdb import RankingDatabase
from pyscenic.genesig import GeneSignature


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

def test_enrichment():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[0]
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    df = enrichment(db, gs)

# -*- coding: utf-8 -*-

import pytest
from pkg_resources import resource_filename

from pyscenic.genesig import GeneSignature
from pyscenic.rnkdb import ParquetRankingDatabase as RankingDatabase

TEST_DATABASE_FNAME = resource_filename('resources.tests', "hg19-tss-centered-10kb-10species.mc9nr.parquet")
TEST_DATABASE_NAME = "hg19-tss-centered-10kb-10species.mc9nr"
TEST_SIGNATURE_FNAME = resource_filename('resources.tests', "c6.all.v6.1.symbols.gmt")

##################################################
# temporarily disable testing of the parqet databases until we have a working test database
from os import path
pytestmark = pytest.mark.skipif(not path.exists(TEST_DATABASE_FNAME), reason="Parquet testing is temporarily disabled.")
##################################################

@pytest.fixture
def db():
    return RankingDatabase(TEST_DATABASE_FNAME, TEST_DATABASE_NAME)


@pytest.fixture
def gs():
    return GeneSignature.from_gmt(TEST_SIGNATURE_FNAME, gene_separator="\t", field_separator="\t", )[0]


def test_init(db):
    assert db.name == TEST_DATABASE_NAME


def test_total_genes(db):
    assert db.total_genes == 22284


def test_genes(db):
    assert len(db.genes) == 22284


def test_load_full(db):
    rankings = db.load_full()
    assert len(rankings.index) == 5
    assert len(rankings.columns) == 22284


def test_load(db, gs):
    rankings = db.load(gs)
    assert len(rankings.index) == 5
    assert len(rankings.columns) == 29

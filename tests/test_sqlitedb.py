# -*- coding: utf-8 -*-

import pytest
from pyscenic.rnkdb import SQLiteRankingDatabase as RankingDatabase
from pyscenic.genesig import GeneSignature
from pkg_resources import resource_filename



NOMENCLATURE = "HGNC"
TEST_DATABASE_FNAME = resource_filename('resources', "hg19-tss-centered-5kb-10species.mc9nr.db")
TEST_DATABASE_NAME = "hg19-tss-centered-5kb-10species"
TEST_SIGNATURE_FNAME = resource_filename('resources', "c6.all.v6.1.symbols.gmt")


@pytest.fixture
def db():
    return RankingDatabase(TEST_DATABASE_FNAME, TEST_DATABASE_NAME, NOMENCLATURE)

@pytest.fixture
def gs():
    return GeneSignature.from_gmt(TEST_SIGNATURE_FNAME, NOMENCLATURE,
                                  gene_separator="\t", field_separator="\t", )[0]


def test_init(db):
    assert db.name == TEST_DATABASE_NAME
    assert db.nomenclature == NOMENCLATURE

def test_total_genes(db):
    assert db.total_genes == 29

def test_load_full(db):
    rankings = db.load_full()
    assert len(rankings.index) == 24453
    assert len(rankings.columns) == 29

def test_load(db, gs):
    rankings = db.load(gs)
    assert len(rankings.index) == 24453
    assert len(rankings.columns) == 29

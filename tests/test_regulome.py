# -*- coding: utf-8 -*-

from pyscenic.regulome import module2regulome

import os, yaml
from configparser import ConfigParser
from pyscenic.rnkdb import SQLiteRankingDatabase as RankingDatabase
from pyscenic.genesig import GeneSignature, Regulome
from pyscenic.utils import load_motif_annotations

TEST_DATABASE = "hg19-500bp-upstream-10species"
TEST_SIGNATURE = "msigdb_cancer_c6"

RESOURCES_FOLDER="/Users/bramvandesande/Projects/lcb/resources"
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")


def load_db_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_sqlitedb.ini'))
    return config[section]

def load_gs_info(section):
    config = ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'test_genesig.ini'))
    return config[section]

def test_module2regulome():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[3]
    db = RankingDatabase(**load_db_info(TEST_DATABASE))
    module = Regulome(gs.name, gs.nomenclature, gs.gene2weights, "TP53")
    motif_annotations = load_motif_annotations(MOTIF_ANNOTATIONS_FNAME)
    reg = module2regulome(db, module, motif_annotations)

def test_to_from_yaml():
    gs = GeneSignature.from_gmt(gene_separator="\t", field_separator="\t", **load_gs_info(TEST_SIGNATURE))[3]
    regulome = Regulome(gs.name, gs.nomenclature, gs.gene2weights, "TP53")
    yml = yaml.dump(regulome)
    regulome2 = yaml.load(yml)
    assert regulome.name == regulome2.name
    assert regulome.nomenclature == regulome2.nomenclature
    assert regulome.genes == regulome2.genes
    assert regulome.weights == regulome2.weights
    assert regulome.score == regulome2.score
    assert regulome.context == regulome2.context


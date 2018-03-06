# -*- coding: utf-8 -*-

from pyscenic.genesig import GeneSignature, Regulome
from configparser import ConfigParser
import os
import pytest
import attr
from pkg_resources import resource_filename


TEST_SIGNATURE = "msigdb_cancer_c6"
TEST_SIGNATURE_FNAME = resource_filename('resources', "c6.all.v6.1.symbols.gmt.txt")
NOMENCLATURE = "HGNC"


def test_init1():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert gs1.nomenclature == 'HGNC'
    assert len(gs1) == 2
    assert gs1.gene2weights['TP53'] == 1.0
    assert gs1.gene2weights['SOX4'] == 1.0


def test_init2():
    gs1 = GeneSignature(name="test1", gene2weights=[('TP53', 0.5), ('SOX4', 0.75)], nomenclature="HGNC")
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert gs1.nomenclature == 'HGNC'
    assert len(gs1) == 2
    assert gs1.gene2weights['TP53'] == 0.5
    assert gs1.gene2weights['SOX4'] == 0.75


def test_init3():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.5, 'SOX4': 0.75}, nomenclature="HGNC")
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert gs1.nomenclature == 'HGNC'
    assert len(gs1) == 2
    assert gs1.gene2weights['TP53'] == 0.5
    assert gs1.gene2weights['SOX4'] == 0.75
    

def test_immut():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.5, 'SOX4': 0.75}, nomenclature="HGNC")
    with pytest.raises(attr.exceptions.FrozenInstanceError):
        gs1.name = 'rename'
    with pytest.raises(attr.exceptions.FrozenInstanceError):
        gs1.nomenclature = 'MGI'
    with pytest.raises(TypeError):
        gs1.gene2weights['TP53'] = 0.6


def test_genes():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.5, 'SOX4': 0.75}, nomenclature="HGNC")
    assert gs1.genes == ('SOX4', 'TP53')


def test_dict():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.5, 'SOX4': 0.75}, nomenclature="HGNC")
    assert gs1['TP53'] == 0.5
    assert gs1['SOX4'] == 0.75


def test_rename():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.5, 'SOX4': 0.75}, nomenclature="HGNC")
    gs2 = gs1.rename('test2')
    assert 'TP53' in gs2
    assert 'SOX4' in gs2
    assert gs2.name == 'test2'
    assert gs2.nomenclature == 'HGNC'
    assert len(gs2) == 2
    assert gs2.gene2weights['TP53'] == 0.5
    assert gs2.gene2weights['SOX4'] == 0.75


def test_union1():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX2'], nomenclature="HGNC")
    gsu = gs1.union(gs2)
    assert 'TP53' in gsu
    assert 'SOX4' in gsu
    assert 'SOX2' in gsu
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 3


def test_union2():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['Tp53', 'Sox4'], nomenclature="MGI")
    with pytest.raises(AssertionError):
        gsu = gs1.union(gs2)


def test_union3():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.8, 'SOX4': 0.75}, nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights={'TP53': 0.3, 'SOX2': 0.60}, nomenclature="HGNC")
    gsu = gs1.union(gs2)
    assert 'TP53' in gsu
    assert gsu.gene2weights['TP53'] == 0.8
    assert 'SOX4' in gsu
    assert gsu.gene2weights['SOX4'] == 0.75
    assert 'SOX2' in gsu
    assert gsu.gene2weights['SOX2'] == 0.6
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 3


def test_diff1():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX2'], nomenclature="HGNC")
    gsu = gs1.difference(gs2)
    assert 'SOX4' in gsu
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 1


def test_diff2():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['Tp53', 'Sox4'], nomenclature="MGI")
    with pytest.raises(AssertionError):
        gsu = gs1.difference(gs2)


def test_diff3():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.8, 'SOX4': 0.75}, nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights={'TP53': 0.3, 'SOX2': 0.60}, nomenclature="HGNC")
    gsu = gs1.difference(gs2)
    assert 'SOX4' in gsu
    assert gsu.gene2weights['SOX4'] == 0.75
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 1


def test_intersection1():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX2'], nomenclature="HGNC")
    gsu = gs1.intersection(gs2)
    assert 'TP53' in gsu
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 1


def test_intersection2():
    gs1 = GeneSignature(name="test1", gene2weights=['TP53', 'SOX4'], nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights=['Tp53', 'Sox4'], nomenclature="MGI")
    with pytest.raises(AssertionError):
        gsu = gs1.intersection(gs2)


def test_intersection3():
    gs1 = GeneSignature(name="test1", gene2weights={'TP53': 0.8, 'SOX4': 0.75}, nomenclature="HGNC")
    gs2 = GeneSignature(name="test1", gene2weights={'TP53': 0.3, 'SOX2': 0.60}, nomenclature="HGNC")
    gsu = gs1.intersection(gs2)
    assert gsu.nomenclature == 'HGNC'
    assert len(gsu) == 1
    assert 'TP53' in gsu
    assert gsu.gene2weights['TP53'] == 0.8


def test_regulome():
    reg = Regulome(name='TP53 regulome', gene2weights={'TP53': 0.8, 'SOX4': 0.75}, nomenclature="HGNC", transcription_factor="TP53")
    assert reg.transcription_factor == "TP53"


def test_load_gmt():
    gss = GeneSignature.from_gmt(field_separator='\t', gene_separator='\t', nomenclature=NOMENCLATURE,
                                 fname=TEST_SIGNATURE_FNAME)
    # http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6
    assert len(gss) == 189
    assert gss[0].name == "GLI1_UP.V1_DN"
    assert "COPZ1" in gss[0]
    assert len(gss[0]) == 29


def test_add():
    gss = GeneSignature.from_gmt(field_separator='\t', gene_separator='\t', nomenclature=NOMENCLATURE,
                                 fname=TEST_SIGNATURE_FNAME)
    res = gss[0].add("MEF2")
    assert "MEF2" in res


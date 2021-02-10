# -*- coding: utf-8 -*-

import attr
import pytest
from pkg_resources import resource_filename

from pyscenic.genesig import GeneSignature, Regulon

TEST_SIGNATURE = "msigdb_cancer_c6"
TEST_SIGNATURE_FNAME = resource_filename('resources.tests', "c6.all.v6.1.symbols.gmt")


def test_init1():
    gs1 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX4'])
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert len(gs1) == 2
    assert gs1.gene2weight['TP53'] == 1.0
    assert gs1.gene2weight['SOX4'] == 1.0


def test_init2():
    gs1 = GeneSignature(name="test1", gene2weight=[('TP53', 0.5), ('SOX4', 0.75)])
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert len(gs1) == 2
    assert gs1.gene2weight['TP53'] == 0.5
    assert gs1.gene2weight['SOX4'] == 0.75


def test_init3():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.5, 'SOX4': 0.75})
    assert 'TP53' in gs1
    assert 'SOX4' in gs1
    assert gs1.name == 'test1'
    assert len(gs1) == 2
    assert gs1.gene2weight['TP53'] == 0.5
    assert gs1.gene2weight['SOX4'] == 0.75


def test_immut():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.5, 'SOX4': 0.75})
    with pytest.raises(attr.exceptions.FrozenInstanceError):
        gs1.name = 'rename'
    with pytest.raises(TypeError):
        gs1.gene2weight['TP53'] = 0.6


def test_genes():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.5, 'SOX4': 0.75})
    assert gs1.genes == ('SOX4', 'TP53')


def test_dict():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.5, 'SOX4': 0.75})
    assert gs1['TP53'] == 0.5
    assert gs1['SOX4'] == 0.75


def test_rename():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.5, 'SOX4': 0.75})
    gs2 = gs1.rename('test2')
    assert 'TP53' in gs2
    assert 'SOX4' in gs2
    assert gs2.name == 'test2'
    assert len(gs2) == 2
    assert gs2.gene2weight['TP53'] == 0.5
    assert gs2.gene2weight['SOX4'] == 0.75


def test_union1():
    gs1 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX4'])
    gs2 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX2'])
    gsu = gs1.union(gs2)
    assert 'TP53' in gsu
    assert 'SOX4' in gsu
    assert 'SOX2' in gsu
    assert len(gsu) == 3


def test_union3():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.8, 'SOX4': 0.75})
    gs2 = GeneSignature(name="test1", gene2weight={'TP53': 0.3, 'SOX2': 0.60})
    gsu = gs1.union(gs2)
    assert 'TP53' in gsu
    assert gsu.gene2weight['TP53'] == 0.8
    assert 'SOX4' in gsu
    assert gsu.gene2weight['SOX4'] == 0.75
    assert 'SOX2' in gsu
    assert gsu.gene2weight['SOX2'] == 0.6
    assert len(gsu) == 3


def test_diff1():
    gs1 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX4'])
    gs2 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX2'])
    gsu = gs1.difference(gs2)
    assert 'SOX4' in gsu
    assert len(gsu) == 1


def test_diff3():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.8, 'SOX4': 0.75})
    gs2 = GeneSignature(name="test1", gene2weight={'TP53': 0.3, 'SOX2': 0.60})
    gsu = gs1.difference(gs2)
    assert 'SOX4' in gsu
    assert gsu.gene2weight['SOX4'] == 0.75
    assert len(gsu) == 1


def test_intersection1():
    gs1 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX4'])
    gs2 = GeneSignature(name="test1", gene2weight=['TP53', 'SOX2'])
    gsu = gs1.intersection(gs2)
    assert 'TP53' in gsu
    assert len(gsu) == 1


def test_intersection3():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.8, 'SOX4': 0.75})
    gs2 = GeneSignature(name="test1", gene2weight={'TP53': 0.3, 'SOX2': 0.60})
    gsu = gs1.intersection(gs2)
    assert len(gsu) == 1
    assert 'TP53' in gsu
    assert gsu.gene2weight['TP53'] == 0.8


def test_regulon():
    reg = Regulon(name='TP53 regulon', gene2weight={'TP53': 0.8, 'SOX4': 0.75}, transcription_factor="TP53", gene2occurrence={"TP53": 1})
    assert reg.transcription_factor == "TP53"


def test_load_gmt():
    gss = GeneSignature.from_gmt(field_separator='\t', gene_separator='\t', fname=TEST_SIGNATURE_FNAME)
    # http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6
    assert len(gss) == 189
    assert gss[0].name == "GLI1_UP.V1_DN"
    assert "COPZ1" in gss[0]
    assert len(gss[0]) == 29


def test_add():
    gss = GeneSignature.from_gmt(field_separator='\t', gene_separator='\t', fname=TEST_SIGNATURE_FNAME)
    res = gss[0].add("MEF2")
    assert "MEF2" in res


def test_noweights():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.8, 'SOX4': 0.75})
    gs2 = gs1.noweights()
    assert gs1['TP53'] == 0.8
    assert gs2['TP53'] == 1.0

    reg1 = Regulon(name='TP53 regulon', gene2weight={'TP53': 0.8, 'SOX4': 0.75}, transcription_factor="TP53", gene2occurrence={"TP53": 1})
    reg2 = reg1.noweights()
    assert reg1['TP53'] == 0.8
    assert reg2['TP53'] == 1.0
    assert isinstance(reg2, Regulon)


def test_head():
    gs1 = GeneSignature(name="test1", gene2weight={'TP53': 0.8, 'SOX4': 0.75})
    gs2 = gs1.head(1)
    assert gs2['TP53'] == 0.8
    assert len(gs2) == 1
    gs2 = gs1.head(2)
    assert gs2['TP53'] == 0.8
    assert gs2['SOX4'] == 0.75
    assert len(gs2) == 2

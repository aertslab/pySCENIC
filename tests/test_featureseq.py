# coding=utf-8


from pyscenic.featureseq import Feature, FeatureSeq


def test_feature_overlap():
    f1 = Feature.from_string('chr1 12 50 feature1 10.0 +')
    f2 = Feature.from_string('chr1 40 60 feature2 10.0 -')
    f3 = Feature.from_string('chr1 60 80 feature2 10.0 -')
    f4 = Feature.from_string('chr2 60 80 feature2 10.0 -')
    assert f1.has_overlap_with(f2)
    assert not f1.has_overlap_with(f3)
    assert not f3.has_overlap_with(f4)
    assert not f2.has_overlap_with(f3)


def test_feature_contains():
    f1 = Feature.from_string('chr1 50 60 feature1 10.0 +')
    f2 = Feature.from_string('chr1 40 60 feature2 10.0 -')
    f3 = Feature.from_string('chr2 40 80 feature2 10.0 -')
    assert f1 in f2
    assert not f1 in f3


def test_feature_bp_overlap():
    f1 = Feature.from_string('chr1 12 50 feature1 10.0 +')
    f2 = Feature.from_string('chr1 40 60 feature2 10.0 -')
    f3 = Feature.from_string('chr1 60 80 feature2 10.0 -')
    f4 = Feature.from_string('chr2 60 80 feature2 10.0 -')
    assert f1.get_overlap_in_bp_with(f2) == 10
    assert f1.get_overlap_in_bp_with(f3) == 0
    assert f2.get_overlap_in_bp_with(f3) == 0
    assert f3.get_overlap_in_bp_with(f4) == 0
    assert f3.get_overlap_in_bp_with(f3) == len(f3)



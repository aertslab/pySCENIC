# coding=utf-8

import gzip
import os
from operator import attrgetter
from functools import partial
from enum import Enum
from pkg_resources import resource_stream
from typing import Tuple, Type
from .rnkdb import InvertedRankingDatabase
from .featureseq import FeatureSeq
from .genesig import GeneSignature
from cytoolz import mapcat


REGIONS_BED_EXTENSION = "bed.gz"
REGION_NOMENCLATURE = "regions"


def _strip_suffix(s):
    return s.split('#')[0]


class Delineation(Enum):
    HG19_500BP_UP = FeatureSeq.from_bed_file(
        gzip.open(
            resource_stream('resources.delineations', "hg19-limited-upstream500.bed.gz")),
            transform=_strip_suffix)
    HG19_5KBP_CTR = FeatureSeq.from_bed_file(
        gzip.open(
            resource_stream('resources.delineations', "hg19-limited-upstream5000-tss-downstream5000-full-transcript.bed.gz")),
            transform=_strip_suffix)
    HG19_10KBP_CTR = FeatureSeq.from_bed_file(
        gzip.open(
            resource_stream('resources.delineations', "hg19-limited-upstream10000-tss-downstream10000-full-transcript.bed.gz")),
        transform=_strip_suffix)


class RegionRankingDatabase(InvertedRankingDatabase):
    def __init__(self, fname: str, name: str, nomenclature: str = REGION_NOMENCLATURE):
        super().__init__(fname, name, nomenclature)

        basename = os.path.basename(fname).split('.')[0]
        dirname = os.path.dirname(fname)
        region_delineation_fname = os.path.join(dirname, '{}.{}'.format(basename, REGIONS_BED_EXTENSION))
        with gzip.open(region_delineation_fname) as f:
            self._regions = FeatureSeq.from_bed_file(f)

    @property
    def regions(self) -> FeatureSeq:
        return self._regions

    @property
    def identifiers(self) -> Tuple[str]:
        return self.genes


def convert(sig: Type[GeneSignature], db: RegionRankingDatabase, delineation: Delineation, fraction: float = 0.80) -> Type[GeneSignature]:
    """
    Convert a signature of gene symbols to a signature of region identifiers.

    :param sig: The signature of gene symbols.
    :param db: The region database.
    :param delineation: The regulatory region delineation for genes.
    :param fraction: The fraction of overlap to take into account.
    :return: The signature of region identifiers.
    """
    assert sig
    assert db
    assert delineation

    region_identifiers = frozenset(
        map(attrgetter('name'),
            mapcat(list,
                map(partial(db.regions.intersection, fraction=fraction),
                   map(delineation.get, sig.genes)))))
    return sig.copy(gene2weight=region_identifiers, nomenclature=REGION_NOMENCLATURE)





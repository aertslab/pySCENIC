# coding=utf-8

import gzip
import os
from operator import attrgetter
from enum import Enum
from pkg_resources import resource_stream
from typing import Tuple, Type
from .rnkdb import InvertedRankingDatabase
from .featureseq import FeatureSeq
from .genesig import GeneSignature
from cytoolz import merge_with, memoize
from itertools import repeat



REGIONS_BED_EXTENSION = "bed.gz"
REGION_NOMENCLATURE = "regions"


# TODO: Dynamic creation of candidate regulatory regions starting from a gene annotation file and a collection of parameters.
# TODO: would greatly improve the additional functionality provided by the region-based approach.
class Delineation(Enum):
    HG19_500BP_UP = "hg19-limited-upstream500.bed.gz"
    HG19_5KBP_CTR = "hg19-limited-upstream5000-tss-downstream5000-full-transcript.bed.gz"
    HG19_10KBP_CTR = "hg19-limited-upstream10000-tss-downstream10000-full-transcript.bed.gz"


@memoize
def load(delineation: Delineation) -> FeatureSeq:
    def strip_suffix(s):
        return s.split('#')[0]

    return FeatureSeq.from_bed_file(
                gzip.open(
                    resource_stream('resources.delineations', delineation.value), "rt"),
                        transform=strip_suffix)


class RegionRankingDatabase(InvertedRankingDatabase):
    def __init__(self, fname: str, name: str, nomenclature: str = REGION_NOMENCLATURE):
        super().__init__(fname, name, nomenclature)

        basename = os.path.basename(fname).split('.')[0]
        dirname = os.path.dirname(fname)
        region_delineation_fname = os.path.join(dirname, '{}.{}'.format(basename, REGIONS_BED_EXTENSION))
        with gzip.open(region_delineation_fname, "rt") as f:
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

    # Every gene is transformed into a dictionary that maps the name of a feature to the weight of the corresponding gene.
    # These mappings are then combined taking the maximum of multiple values exists for a key.
    identifier2weight = merge_with(max, (dict(zip(
                                            map(attrgetter('name'), db.regions.intersection(load(delineation).get(gene),
                                                                                            fraction=fraction)),
                                            repeat(weight)))
                                         for gene, weight in sig.gene2weight.items()))

    #region_identifiers = frozenset(
    #    map(attrgetter('name'),
    #        mapcat(list,
    #            map(partial(db.regions.intersection, fraction=fraction),
    #               map(load(delineation).get, sig.genes)))))
    return sig.copy(gene2weight=identifier2weight, nomenclature=REGION_NOMENCLATURE)





# coding=utf-8

import gzip
import os
from operator import attrgetter
from enum import Enum
from pkg_resources import resource_stream
from typing import Tuple, Type, Iterator, Sequence
from .rnkdb import InvertedRankingDatabase
from .featureseq import FeatureSeq
from .genesig import GeneSignature, Regulon
from cytoolz import merge_with, memoize
from itertools import repeat, chain
import pandas as pd

REGIONS_BED_EXTENSION = "bed.gz"


# TODO: Dynamic creation of candidate regulatory regions starting from a gene annotation file and a collection of parameters.
# would greatly improve the additional functionality provided by the region-based approach.
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
    # TODO: Improve read performance of inverted design
    # the inverted database design can significantly reduce the size on disk (from 120Gb to 4,7Gb for the 1M regions-24K
    # features human database) but with a huge impact on read performance (from 1 second to several minutes for a typical
    # signature). Potential mitigation challenges are:
    #  - Recoding the decompression of the inverted design using C++ and its standard template library. Impact will only
    #    be a constant factor while personal investment will be substantial.
    #  - Break the clean interface between database storage and AUC and recovery curve calculation. The name of the
    #    genes that are at the top of a feature ranking is not necessary to calculate the AUC for enrichment of this
    #    feature. This fact can be used so that decompression times could be significantly reduced.
    #
    def __init__(self, fname: str, name: str):
        super().__init__(fname, name)

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

    return sig.copy(gene2weight=identifier2weight,
                    context=frozenset(chain(sig.context, [str(delineation)])))


def df2regulons(df: pd.DataFrame, gene_regulons: Sequence[Regulon],
                db: RegionRankingDatabase,
                delineations: Iterator[Delineation] = Delineation) -> Sequence[Regulon]:
    """
    Create gene-based regulons from a dataframe of enriched motifs.
    """

    # TODO: This method needs to be implemented to have a fully functional regions-based pipeline.
    # There are many ways to provide the link from enhancer-based regulons back to the original and pruned gene-based
    # regulons. From a high level there are two approaches: either the software does the bookkeeping (i.e. we keep
    # the link between regions and their corresponding genes when converting from gene to region-based regulons) or
    # the end user has to this himselves. The current function signature implements the latter strategy, i.e. the
    # user has to supply the candidate regulatory regions for genes and the original gene-based regulons to
    # be able to do the conversion.

    # Algorithm's skeleton: For each row:
    # 1. Find the appropriate delineation and convert regions in the TargetGenes column of the dataframe to their
    #    corresponding geneIDs.
    # 1.1. This mapping between a region ID and (multiple) gene IDs needs the regions property of a RegionRankingDatabase
    #      together with a regulatory region delineation which can be derived from an element in the context field.
    #      Caveat: (1) Either we provide the delineation instance itself or a string representation. Currently the later
    #      is being done and therefore we need the list of all potential delineations as an extra input parameter.
    #      (2) It is better the precalculate or at least cache this operation because it will be needed for many more
    #      rows.
    # 2. Filter the list of gene IDs based on presence in the original corresponding gene-based regulon.
    # 2.1 The main challenge is to find the corresponding gene-based regulon for each row in the enriched motif dataframe.
    #     This mapping can be resolved by using the name of transcription factor together with elements provided in the
    #     context column (mainly method used for defining a module from the adjacencies; the name of the database is
    #     less informational).

    raise NotImplemented


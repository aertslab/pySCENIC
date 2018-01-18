# -*- coding: utf-8 -*-

import attr
import re
from typing import Mapping, List
from frozendict import frozendict
from collections.abc import Iterable, Mapping
from itertools import repeat
from cytoolz import merge_with, dissoc, keyfilter, first, second


def convert(genes):
    # Genes supplied as dictionary.
    if isinstance(genes, Mapping):
        return frozendict(genes)
    # Genes supplied as iterable of (gene, weight) tuples.
    elif isinstance(genes, Iterable) and all(isinstance(n, tuple) for n in genes):
        return frozendict(genes)
    # Genes supplied as iterable of genes.
    elif isinstance(genes, Iterable) and all(isinstance(n, str) for n in genes):
        return frozendict(zip(genes, repeat(1.0)))


@attr.s(frozen=True)
class GeneSignature:
    """
    A class of gene signatures, i.e. a set of genes.
    """

    @classmethod
    def from_gmt(cls, fname: str, nomenclature: str, field_separator: str =';', gene_separator=';') -> List['GeneSignature']:
        """

        :param fname:
        :param field_separator:
        :param gene_separator:
        :return:
        """
        # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

        def signatures():
            with open(fname, "r") as file:
                for line in file:
                    if line.startswith("#") or not line.strip():
                        continue
                    columns = re.split(field_separator, line.rstrip())
                    genes = columns[2:] if field_separator == gene_separator else columns[2].split(gene_separator)
                    yield GeneSignature(name=columns[0], nomenclature=nomenclature, genes=genes)
        return list(signatures())

    @classmethod
    def from_grp(cls, fname) -> 'GeneSignature':
        """

        :param fname:
        :return:
        """

        # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
        #TODO
        pass

    name: str = attr.ib()
    nomenclature: str = attr.ib(default="HGNC")
    gene2weights: Mapping[str, float] = attr.ib(converter=convert)

    @name.validator
    def name_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have a non-empty name.")

    @nomenclature.validator
    def name_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have a nomenclature.")

    @gene2weights.validator
    def gene2weights_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have at least one gene.")

    @property
    def genes(self):
        """
        Return genes in this signature. Genes are sorted in descending order according to weight.
        """
        return tuple(map(first, sorted(self.gene2weights.items(), key=second, reverse=True)))

    def rename(self, name: str) -> 'GeneSignature':
        """
        Rename this signature.

        :param name: The new name.
        :return: the new :class:`GeneSignature` instance.
        """
        return GeneSignature(name=name, nomenclature=self.nomenclature, genes=self.genes)

    def union(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the union of this signature and the other supplied
        signature.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Union of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} | {})".format(self.name, other.name),
                             nomenclature=self.nomenclature,
                             genes=frozendict(merge_with(max, self.gene2weights, other.gene2weights)))

    def difference(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the difference of this signature and the supplied other
        signature.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Difference of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} - {})".format(self.name, other.name),
                         nomenclature=self.nomenclature,
                         genes=frozendict(dissoc(self.gene2weights, other.gene2weights.keys())))

    def intersection(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the intersection of this signature and the supplied other
        signature.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Intersection of gene signatures is only possible when both signatures use same nomenclature for genes."
        genes = set(self.gene2weights.keys()).intersection(set(other.gene2weights.keys()))
        return GeneSignature(name="({} & {})".format(self.name, other.name),
                             nomenclature=self.nomenclature,
                             genes=frozenset(keyfilter(lambda k: k in genes,
                                                       merge_with(max, self.gene2weights, other.gene2weights))))

    def __len__(self):
        """
        The number of genes in this signature.
        """
        return len(self.genes)


@attr.s(frozen=True)
class Regulome(GeneSignature):
    """
    A regulome is a gene signature that defines the target genes of a transcription factor.
    """
    transcription_factor: str = attr.ib()

    @transcription_factor.validator
    def non_empty(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A regulome must have a transcription factor.")

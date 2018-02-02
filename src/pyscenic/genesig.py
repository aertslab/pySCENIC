# -*- coding: utf-8 -*-

import re
import os
from collections.abc import Iterable, Mapping
from itertools import repeat
from typing import Mapping, List, Tuple

import attr
from cytoolz import merge_with, dissoc, keyfilter, first, second
from frozendict import frozendict

from cytoolz import memoize


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
    def from_gmt(cls, fname: str, nomenclature: str, field_separator: str =',', gene_separator=',') -> List['GeneSignature']:
        """
        Load gene signatures from a GMT file.

        :param fname: The filename.
        :param nomenclature: The nomenclature of the genes.
        :param field_separator: The separator that separates fields in a line.
        :param gene_separator: The separator that separates the genes.
        :return: A list of signatures.
        """
        # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
        assert os.path.exists(fname), "{} does not exist.".format(fname)

        def signatures():
            with open(fname, "r") as file:
                for line in file:
                    if line.startswith("#") or not line.strip():
                        continue
                    columns = re.split(field_separator, line.rstrip())
                    genes = columns[2:] if field_separator == gene_separator else columns[2].split(gene_separator)
                    yield GeneSignature(name=columns[0], nomenclature=nomenclature, gene2weights=genes)
        return list(signatures())

    @classmethod
    def from_grp(cls, fname, name: str, nomenclature: str) -> 'GeneSignature':
        """
        Load gene signature from GRP file.

        :param fname: The filename.
        :param name: The name of the resulting signature.
        :param nomenclature: The nomenclature of the genes.
        :return: A signature.
        """
        # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
        assert os.path.exists(fname), "{} does not exist.".format(fname)
        with open(fname, "r") as file:
            return GeneSignature(name=name,
                             nomenclature=nomenclature,
                             gene2weights=[line.rstrip() for line in file if not line.startswith("#") and line.strip()])

    @classmethod
    def from_rnk(cls, fname: str, name: str, nomenclature: str, field_separator=",") -> 'GeneSignature':
        """
        Reads in a signature from an RNK file. This format associates weights with the genes part of the signature.

        :param fname: The filename.
        :param name: The name of the resulting signature.
        :param nomenclature: The nomenclature of the genes.
        :param field_separator: The separator that separates fields in a line.
        :return: A signature.
        """
        # https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
        assert os.path.exists(fname), "{} does not exist.".format(fname)

        def columns():
            with open(fname, "r") as file:
                for line in file:
                    if line.startswith("#") or not line.strip():
                        continue
                    columns = tuple(map(str.rstrip, re.split(field_separator, line)))
                    assert len(columns) == 2, "Invalid file format."
                    yield columns

        return GeneSignature(name=name,
                             nomenclature=nomenclature,
                             gene2weights=list(columns()))

    name: str = attr.ib()
    nomenclature: str = attr.ib()
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
    @memoize
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
        return GeneSignature(name=name, nomenclature=self.nomenclature, gene2weights=self.gene2weights)

    def union(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the union of this signature and the other supplied
        signature.

        The weight associated with the genes in the intersection is the maximum of the weights in the composing signatures.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Union of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} | {})".format(self.name, other.name),
                             nomenclature=self.nomenclature,
                             gene2weights=frozendict(merge_with(max, self.gene2weights, other.gene2weights)))

    def difference(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the difference of this signature and the supplied other
        signature.

        The weight associated with the genes in the difference are taken from this gene signature.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Difference of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} - {})".format(self.name, other.name),
                         nomenclature=self.nomenclature,
                         gene2weights=frozendict(dissoc(dict(self.gene2weights), *other.genes)))

    def intersection(self, other: 'GeneSignature') -> 'GeneSignature':
        """
        Creates a new :class:`GeneSignature` instance which is the intersection of this signature and the supplied other
        signature.

        The weight associated with the genes in the intersection is the maximum of the weights in the composing signatures.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Intersection of gene signatures is only possible when both signatures use same nomenclature for genes."
        genes = set(self.gene2weights.keys()).intersection(set(other.gene2weights.keys()))
        return GeneSignature(name="({} & {})".format(self.name, other.name),
                             nomenclature=self.nomenclature,
                             gene2weights=frozendict(keyfilter(lambda k: k in genes,
                                                       merge_with(max, self.gene2weights, other.gene2weights))))

    def __len__(self):
        """
        The number of genes in this signature.
        """
        return len(self.genes)

    def __contains__(self, item):
        """
        Checks if a gene is part of this signature.
        """
        return item in self.gene2weights.keys()

    def __getitem__(self, item):
        """
        Return the weight associated with a gene.
        """
        return self.gene2weights[item]

    def __str__(self):
        """
        Returns a readable string representation.
        """
        return "[]".format(",".join(self.genes))

    def __repr__(self):
        """
        Returns a unambiguous string representation.
        """
        return "{}(name=\"{}\",nomenclature={},n={},genes=[{}])".format(
            self.__class__.__name__,
            self.name,
            self.nomenclature,
            len(self.genes),
            ",".join(map("\"{}\"".format,self.genes)))


@attr.s(frozen=True)
class Regulome(GeneSignature):
    """
    A regulome is a gene signature that defines the target genes of a transcription factor.
    """
    transcription_factor: str = attr.ib()
    context: Tuple[str] = attr.ib(default=())
    score: float = attr.ib(default=0.0)

    @transcription_factor.validator
    def non_empty(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A regulome must have a transcription factor.")

# -*- coding: utf-8 -*-

import re
import os
from collections.abc import Iterable, Mapping
from itertools import repeat
from typing import Mapping, List, FrozenSet, Type

import attr
import yaml
from cytoolz import merge_with, dissoc, keyfilter, first, second
from frozendict import frozendict
from itertools import chain

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
class GeneSignature(yaml.YAMLObject):
    """
    A class of gene signatures, i.e. a set of genes that are biologically related.
    """

    yaml_tag = u'!GeneSignature'

    @classmethod
    def to_yaml(cls, dumper, data):
        dict_representation = {
            'name': data.name,
            'nomenclature': data.nomenclature,
            'genes': list(data.genes),
            'weights': list(data.weights)
        }
        return dumper.represent_mapping(cls.yaml_tag,
                                        dict_representation,
                                            cls)

    @classmethod
    def from_yaml(cls, loader, node):
        data = loader.construct_mapping(node, cls)
        return GeneSignature(name=data['name'],
                             nomenclature=data['nomenclature'],
                             gene2weight=zip(data['genes'], data['weights']))

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
                    yield GeneSignature(name=columns[0], nomenclature=nomenclature, gene2weight=genes)
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
                             gene2weight=[line.rstrip() for line in file if not line.startswith("#") and line.strip()])

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
                             gene2weight=list(columns()))

    name = attr.ib()  # str
    nomenclature = attr.ib()  # str
    gene2weight = attr.ib(converter=convert)  # Mapping[str, float]

    @name.validator
    def name_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have a non-empty name.")

    @nomenclature.validator
    def name_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have a nomenclature.")

    @gene2weight.validator
    def gene2weight_validator(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A gene signature must have at least one gene.")

    @property
    @memoize
    def genes(self):
        """
        Return genes in this signature. Genes are sorted in descending order according to weight.
        """
        return tuple(map(first, sorted(self.gene2weight.items(), key=second, reverse=True)))

    @property
    @memoize
    def weights(self):
        """
        Return the weights of the genes in this signature. Genes are sorted in descending order according to weight.
        """
        return tuple(map(second, sorted(self.gene2weight.items(), key=second, reverse=True)))

    def rename(self, name: str) -> Type['GeneSignature']:
        """
        Rename this signature.

        :param name: The new name.
        :return: the new :class:`GeneSignature` instance.
        """
        return GeneSignature(name=name, nomenclature=self.nomenclature, gene2weight=self.gene2weight)

    def add(self, gene_symbol, weight=1.0)  -> Type['GeneSignature']:
        """
        Add an extra gene symbol to this signature.
        :param gene_symbol: The symbol of the gene.
        :param weight: The weight.
        :return: the new :class:`GeneSignature` instance.
        """
        return GeneSignature(name=self.name, nomenclature=self.nomenclature,
                             gene2weight=list(chain(self.gene2weight.items(), [(gene_symbol, weight)])))

    def _union_impl(self, other):
        return frozendict(merge_with(max, self.gene2weight, other.gene2weight))

    def union(self, other: Type['GeneSignature']) -> Type['GeneSignature']:
        """
        Creates a new :class:`GeneSignature` instance which is the union of this signature and the other supplied
        signature.

        The weight associated with the genes in the intersection is the maximum of the weights in the composing signatures.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Union of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} | {})".format(self.name, other.name) if self.name != other.name else self.name,
                             nomenclature=self.nomenclature,
                             gene2weight=self._union_impl(other))

    def _difference_impl(self, other):
        return frozendict(dissoc(dict(self.gene2weight), *other.genes))

    def difference(self, other: Type['GeneSignature']) -> Type['GeneSignature']:
        """
        Creates a new :class:`GeneSignature` instance which is the difference of this signature and the supplied other
        signature.

        The weight associated with the genes in the difference are taken from this gene signature.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Difference of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} - {})".format(self.name, other.name) if self.name != other.name else self.name,
                         nomenclature=self.nomenclature,
                         gene2weight=self._difference_impl(other))

    def _intersection_impl(self, other):
        genes = set(self.gene2weight.keys()).intersection(set(other.gene2weight.keys()))
        return frozendict(keyfilter(lambda k: k in genes,
                                    merge_with(max, self.gene2weight, other.gene2weight)))

    def intersection(self, other: Type['GeneSignature']) -> Type['GeneSignature']:
        """
        Creates a new :class:`GeneSignature` instance which is the intersection of this signature and the supplied other
        signature.

        The weight associated with the genes in the intersection is the maximum of the weights in the composing signatures.

        :param other: The other :class:`GeneSignature`.
        :return: the new :class:`GeneSignature` instance.
        """
        assert self.nomenclature == other.nomenclature, "Intersection of gene signatures is only possible when both signatures use same nomenclature for genes."
        return GeneSignature(name="({} & {})".format(self.name, other.name) if self.name != other.name else self.name,
                             nomenclature=self.nomenclature,
                             gene2weight=self._intersection_impl(other))

    def noweights(self):
        """
        Create a new gene signature with uniform weights, i.e. all weights are equal and set to 1.0.
        """
        return GeneSignature(name=self.name, nomenclature=self.nomenclature,
                             gene2weight=self.genes)

    def __len__(self):
        """
        The number of genes in this signature.
        """
        return len(self.genes)

    def __contains__(self, item):
        """
        Checks if a gene is part of this signature.
        """
        return item in self.gene2weight.keys()

    def __getitem__(self, item):
        """
        Return the weight associated with a gene.
        """
        return self.gene2weight[item]

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
class Regulon(GeneSignature, yaml.YAMLObject):
    """
    A regulon is a gene signature that defines the target genes of a Transcription Factor (TF) and thereby defines
    a subnetwork of a larger Gene Regulatory Network (GRN) connecting a TF with its target genes.
    """

    yaml_tag = u'!Regulon'

    @classmethod
    def to_yaml(cls, dumper, data):
        dict_representation = {
            'name': data.name,
            'nomenclature': data.nomenclature,
            'genes': list(data.genes),
            'weights': list(data.weights),
            'score': data.score,
            'context': list(data.context),
            'transcription_factor': data.transcription_factor
        }
        return dumper.represent_mapping(cls.yaml_tag,
                                            dict_representation,
                                            cls)

    @classmethod
    def from_yaml(cls, loader, node):
        data = loader.construct_mapping(node, cls)
        return Regulon(name=data['name'],
                         nomenclature=data['nomenclature'],
                         gene2weight=list(zip(data['genes'], data['weights'])),
                         score=data['score'],
                         context=frozenset(data['context']),
                         transcription_factor=data['transcription_factor'])

    transcription_factor = attr.ib()  # str
    context = attr.ib(default=frozenset())  # FrozenSet[str]
    score = attr.ib(default=0.0)  # float

    @transcription_factor.validator
    def non_empty(self, attribute, value):
        if len(value) == 0:
            raise ValueError("A regulon must have a transcription factor.")

    def add(self, gene_symbol, weight=1.0) -> 'Regulon':
        """
        Add an extra gene symbol to this signature.
        :param gene_symbol: The symbol of the gene.
        :param weight: The weight.
        :return: the new :class:`GeneSignature` instance.
        """
        return Regulon(name=self.name, nomenclature=self.nomenclature,
                        context=self.context, transcription_factor=self.transcription_factor,
                        score=self.score,
                             gene2weight=list(chain(self.gene2weight.items(), [(gene_symbol, weight)])))

    def union(self, other: Type['GeneSignature']) -> 'Regulon':
        assert self.nomenclature == other.nomenclature, "Union of gene signatures is only possible when both signatures use same nomenclature for genes."
        assert self.transcription_factor == getattr(other, 'transcription_factor', self.transcription_factor), "Union of two regulons is only possible when same factor."
        return Regulon(name="({} | {})".format(self.name, other.name) if self.name != other.name else self.name,
                             nomenclature=self.nomenclature,
                             transcription_factor=self.transcription_factor,
                             context=self.context.union(getattr(other, 'context', frozenset())),
                             score=max(self.score, getattr(other, 'score', 0.0)),
                             gene2weight=self._union_impl(other))


    def difference(self, other: Type['GeneSignature']) -> 'Regulon':
        assert self.nomenclature == other.nomenclature, "Difference of gene signatures is only possible when both signatures use same nomenclature for genes."
        assert self.transcription_factor == getattr(other, 'transcription_factor', self.transcription_factor), "Difference of two regulons is only possible when same factor."
        return Regulon(name="({} - {})".format(self.name, other.name) if self.name != other.name else self.name,
                             nomenclature=self.nomenclature,
                             transcription_factor=self.transcription_factor,
                             context=self.context.union(getattr(other, 'context', frozenset())),
                             score=max(self.score, getattr(other, 'score', 0.0)),
                             gene2weight=self._difference_impl(other))

    def intersection(self, other: Type['GeneSignature']) -> 'Regulon':
        assert self.nomenclature == other.nomenclature, "Intersection of gene signatures is only possible when both signatures use same nomenclature for genes."
        assert self.transcription_factor == getattr(other, 'transcription_factor', self.transcription_factor), "Intersection of two regulons is only possible when same factor."
        return Regulon(name="({} & {})".format(self.name, other.name) if self.name != other.name else self.name,
                         nomenclature=self.nomenclature,
                         transcription_factor=self.transcription_factor,
                         context=self.context.union(getattr(other, 'context', frozenset())),
                         score=max(self.score, getattr(other, 'score', 0.0)),
                         gene2weight=self._intersection_impl(other))

    def noweights(self):
        """
        Create a new regulon with uniform weights, i.e. all weights are equal and set to 1.0.
        """
        return Regulon(name=self.name,
                        nomenclature=self.nomenclature,
                        transcription_factor=self.transcription_factor,
                        context=self.context,
                        score=self.score,
                         gene2weight=self.genes)

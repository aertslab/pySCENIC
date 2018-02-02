# -*- coding: utf-8 -*-

import abc
import os
import pandas as pd
from .genesig import GeneSignature
from typing import Tuple, Set, Type
from cytoolz import memoize
from abc import ABCMeta, abstractmethod


class RankingDatabase(metaclass=ABCMeta):
    """
    A class of a database of whole genome rankings. The whole genome is ranked for regulatory features of interest, e.g.
    motifs for a transcription factor.
    """

    @classmethod
    def create(cls, fname: str, name: str, nomenclature: str) -> Type['RankingDatabase']:
        """
        Create a ranking database.

        :param fname: The filename of the database.
        :param name: The name of the database.
        :param nomenclature: The nomenclature used for the genes that are being ranked.
        :return: A ranking database.
        """
        assert os.path.isfile(fname), "{} does not exist.".format(fname)
        assert name, "A database should be given a proper name."
        assert nomenclature, "Nomenclature for the genes in a database should be given."

        extension = os.path.splitext(fname)[1]
        from .featherdb import FeatherRankingDatabase
        from .sqlitedb import SQLiteRankingDatabase
        if extension == "feather":
            # noinspection PyTypeChecker
            return FeatherRankingDatabase(fname, name=name, nomenclature=nomenclature)
        elif extension in ("db", "sqlite", "sqlite3"):
            # noinspection PyTypeChecker
            return SQLiteRankingDatabase(fname, name=name, nomenclature=nomenclature)
        else:
            raise ValueError("{} is an unknown extension.".format(extension))

    def __init__(self, fname: str, name: str, nomenclature: str):
        """
        Create a new instance.

        :param fname: The name of the database file.
        :param nomenclature: The gene nomenclature.
        :param name: The name of the database.
        """
        assert os.path.isfile(fname), "Database {0:s} doesn't exist.".format(fname)
        assert name, "Name must be specified."
        assert nomenclature, "Nomenclature must be specified."

        self._fname = fname
        self._name = name
        self._nomenclature = nomenclature

    @property
    def name(self) -> str:
        """
        The name of this database of rankings.
        """
        return self._name

    @property
    def nomenclature(self) -> str:
        """
        The nomenclature used for specifying the genes.
        """
        return self._nomenclature

    @property
    @abstractmethod
    def total_genes(self) -> int:
        """
        The total number of genes ranked.
        """
        pass

    @property
    @abstractmethod
    def genes(self) -> Tuple[str]:
        """
        List of genes ranked according to the regulatory features in this database.
        """
        pass

    @property
    @memoize
    def geneset(self) -> Set[str]:
        """
        Set of genes ranked according to the regulatory features in this database.
        """
        return set(self.genes)

    @abstractmethod
    def load_full(self) -> pd.DataFrame:
        """
        Load the whole database into memory.

        :return: a dataframe.
        """
        pass

    @abstractmethod
    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        """
        Load the ranking of the genes in the supplied signature for all features in this database.

        :param gs: The gene signature.
        :return: a dataframe.
        """
        pass

    def __str__(self):
        """
        Returns a readable string representation.
        """
        return self.name

    def __repr__(self):
        """
        Returns a unambiguous string representation.
        """
        return "{}(name=\"{}\",nomenclature={})".format(
            self.__class__.__name__,
            self.name,
            self.nomenclature)

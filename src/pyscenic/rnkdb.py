# -*- coding: utf-8 -*-

import os
import pandas as pd
from typing import Tuple, Set, Type
from abc import ABCMeta, abstractmethod
import sqlite3
from operator import itemgetter
import numpy as np
from .genesig import GeneSignature
from cytoolz import memoize
from feather.api import write_dataframe, FeatherReader


class RankingDatabase(metaclass=ABCMeta):
    """
    A class of a database of whole genome rankings. The whole genome is ranked for regulatory features of interest, e.g.
    motifs for a transcription factor.

    The rankings of the genes are 0-based.
    """

    def __init__(self, name: str, nomenclature: str):
        """
        Create a new instance.

        :param nomenclature: The gene nomenclature.
        :param name: The name of the database.
        """
        assert name, "Name must be specified."
        assert nomenclature, "Nomenclature must be specified."

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


# SQL query to get the total number of genes in the database.
GENE_ID_COUNT_QUERY = r"SELECT COUNT(*) FROM rankings;"
# SQL query for retrieving the rankings for a particular set of genes.
RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID;"
# SQL query that retrieves the ordered list of features in the database.
FEATURE_IDS_QUERY = r"SELECT motifName FROM motifs ORDER BY idx;"
# SQL query for retrieving the full list of genes scored in this database.
ALL_GENE_IDS_QUERY = r"SELECT geneID FROM rankings ORDER BY geneID;"
# SQL query for retrieving the the whole database.
ALL_RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings ORDER BY geneID;"


class SQLiteRankingDatabase(RankingDatabase):
    """
    A class of a database of whole genome rankings. The whole genome is ranked for regulatory features of interest, e.g.
    motifs for a transcription factor.
    """

    def __init__(self, fname: str, name: str, nomenclature: str):
        """
        Create a new instance.

        :param fname: The name of the SQLite database file.
        :param nomenclature: The gene nomenclature.
        :param name: The name of the database.
        """
        assert os.path.isfile(fname), "Database {0:s} doesn't exist.".format(fname)

        super().__init__(name, nomenclature)

        self._fname = fname
        # Read-only view on SQLite database.
        self._uri = 'file:{}?mode=ro'.format(os.path.abspath(fname))

        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            count = cursor.execute(GENE_ID_COUNT_QUERY).fetchone()
            cursor.close()
        self._gene_count = count[0]

        # Because of problems on same architectures use of unsigned integers is avoided.
        def derive_dtype(n):
            """ Derive datatype for storing 0-based rankings for a given set length. """
            if n <= 2**15:
                # Range int16: -2^15 (= -32768) to 2^15 - 1 (= 32767).
                return np.int16
            else:
                # Range int32: -2^31 (= -2147483648) to 2^31 - 1 (= 2147483647).
                return np.int32
        self._dtype = derive_dtype(self._gene_count)

    @property
    def total_genes(self) -> int:
        """
        The total number of genes ranked.
        """
        return self._gene_count

    @property
    @memoize
    def features(self) -> Tuple[str]:
        """
        List of regulatory features for which whole genome rankings are available in this database.
        """
        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            features = tuple(map(itemgetter(0), cursor.execute(FEATURE_IDS_QUERY).fetchall()))
            cursor.close()
        return features

    @property
    @memoize
    def genes(self) -> Tuple[str]:
        """
        List of genes ranked according to the regulatory features in this database.
        """
        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            genes = tuple(map(itemgetter(0), cursor.execute(ALL_GENE_IDS_QUERY).fetchall()))
            cursor.close()
        return genes

    def load_full(self) -> pd.DataFrame:
        """
        Load the whole database into memory.

        :return: a dataframe.
        """
        # Pre-allocate the matrix.
        rankings = np.empty(shape=(len(self.features), len(self.genes)), dtype=self._dtype)
        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            for idx, (_, ranking) in enumerate(cursor.execute(ALL_RANKINGS_QUERY)):
                rankings[:, idx] = np.frombuffer(ranking, dtype=self._dtype)
            cursor.close()

        return pd.DataFrame(index=self.features, columns=self.genes, data=rankings)

    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        """
        Load the ranking of the genes in the supplied signature for all features in this database.

        :param gs: The gene signature.
        :return: A dataframe.
        """
        assert gs, "A gene signature must be supplied"

        def quoted_csv(values):
            # Escape single quotes (') by using (''), because sometimes ID's contain a single quote.
            def quote(value):
                return "'" + value.replace("'", "''") + "'"
            return ','.join(map(quote, values))

        # For some genes in the signature there might not be a rank available in the database.
        gene_set = self.geneset.intersection(set(gs.genes))
        # Pre-allocate the matrix.
        rankings = np.empty(shape=(len(self.features), len(gene_set)), dtype=self._dtype)
        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            genes = []
            for idx, (gene, ranking) in enumerate(cursor.execute(RANKINGS_QUERY.format(quoted_csv(gene_set)))):
                rankings[:, idx] = np.frombuffer(ranking, dtype=self._dtype)
                genes.append(gene)
            cursor.close()

        return pd.DataFrame(index=self.features, columns=genes, data=rankings)


INDEX_NAME = "features"


class FeatherRankingDatabase(RankingDatabase):
    def __init__(self, fname: str, name: str = None, nomenclature: str = None):
        """
        Create a new feather database.

        :param fname: The filename of the database.
        :param name: The name of the database.
        :param nomenclature: The nomenclature used for the genes that are being ranked.
        """
        assert os.path.isfile(fname), "Database {0:s} doesn't exist.".format(fname)
        super().__init__(name=name, nomenclature=nomenclature)
        # FeatherReader cannot be pickle (important for dask framework) so filename is field instead.
        self._fname = fname

    @property
    @memoize
    def total_genes(self) -> int:
        return FeatherReader(self._fname).num_columns - 1

    @property
    @memoize
    def genes(self) -> Tuple[str]:
        # noinspection PyTypeChecker
        reader = FeatherReader(self._fname)
        return tuple(reader.get_column_name(idx) for idx in range(self.total_genes) if reader.get_column_name(idx) != INDEX_NAME)

    def load_full(self) -> pd.DataFrame:
        return FeatherReader(self._fname).read().set_index(INDEX_NAME)

    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        return FeatherReader(self._fname).read(columns=(INDEX_NAME,) + gs.genes).set_index(INDEX_NAME)


class MemoryDecorator(RankingDatabase):
    """
    A decorator for a ranking database which loads the entire database in memory.
    """
    def __init__(self, db: Type[RankingDatabase]):
        assert db, "Database should be supplied."
        self._db = db
        self._df = db.load_full()
        super().__init__(db._fname, db.name, db.nomenclature)

    @property
    def total_genes(self) -> int:
        return self._db.total_genes

    @property
    def genes(self) -> Tuple[str]:
        return self._db.genes

    def load_full(self) -> pd.DataFrame:
        return self._df

    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        return self._df.loc[:, self._df.columns.isin(gs.genes)]


class DataFrameRankingDatabase(RankingDatabase):
    """
    A ranking database from a dataframe.
    """
    def __init__(self, df: pd.DataFrame, name: str, nomenclature: str):
        self._df = df
        super().__init__(name, nomenclature)

    @property
    def total_genes(self) -> int:
        return len(self._df.columns)

    @property
    def genes(self) -> Tuple[str]:
        return tuple(self._df.columns)

    def load_full(self) -> pd.DataFrame:
        return self._df

    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        return self._df.loc[:, self._df.columns.isin(gs.genes)]

    def save(self, fname: str):
        """
        Save database as feather file.

        :param fname: The name of the file to create.
        """
        assert not os.path.exists(fname)
        df = self._df.copy()
        df.index.name = INDEX_NAME
        df.reset_index(inplace=True) # Index is not stored in feather format. https://github.com/wesm/feather/issues/200
        write_dataframe(df, fname)


def convert2feather(fname: str, out_folder: str, name: str, extension: str="feather") -> str:
    """
    Convert a whole genome rankings database to a feather format based database.

    More information on this format can be found here:
    .. feather-format: https://blog.rstudio.com/2016/03/29/feather/

    :param fname: The filename of the legacy
    :param out_folder: The name of the folder to write the new database to.
    :param name: The name of the rankings database.
    :param extension: The extension of the new database file.
    :return: The filename of the new database.
    """
    assert os.path.isfile(fname), "{} does not exist.".format(fname)
    assert os.path.isdir(out_folder), "{} is not a directory.".format(out_folder)

    feather_fname = os.path.join(out_folder, "{}.{}".format(os.path.splitext(os.path.basename(fname))[0], extension))
    assert not os.path.exists(feather_fname), "{} already exists.".format(feather_fname)

    # Load original database into memory.
    # Caveat: the original storage format of whole genome rankings does not store the metadata, i.e. name and gene
    # nomenclature.
    # The avoid having to specify nomenclature it is set as unknown.
    db = SQLiteRankingDatabase(fname=fname, name=name, nomenclature="UNKNOWN")
    df = db.load_full()
    df.index.name = INDEX_NAME
    df.reset_index(inplace=True) # Index is not stored in feather format. https://github.com/wesm/feather/issues/200
    write_dataframe(df, feather_fname)
    return feather_fname


def open(fname: str, name: str, nomenclature: str) -> Type['RankingDatabase']:
    """
    Open a ranking database.

    :param fname: The filename of the database.
    :param name: The name of the database.
    :param nomenclature: The nomenclature used for the genes that are being ranked.
    :return: A ranking database.
    """
    assert os.path.isfile(fname), "{} does not exist.".format(fname)
    assert name, "A database should be given a proper name."
    assert nomenclature, "Nomenclature for the genes in a database should be given."

    extension = os.path.splitext(fname)[1]
    if extension == "feather":
        # noinspection PyTypeChecker
        return FeatherRankingDatabase(fname, name=name, nomenclature=nomenclature)
    elif extension in ("db", "sqlite", "sqlite3"):
        # noinspection PyTypeChecker
        return SQLiteRankingDatabase(fname, name=name, nomenclature=nomenclature)
    else:
        raise ValueError("{} is an unknown extension.".format(extension))

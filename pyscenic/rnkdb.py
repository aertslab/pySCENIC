# -*- coding: utf-8 -*-
import sqlite3
import os
from operator import itemgetter
import numpy as np
from typing import Tuple
from .genesig import GeneSignature

# SQL query to get the total number of genes in the database.
GENE_ID_COUNT_QUERY = r"SELECT COUNT(*) FROM rankings;"
# SQL query for retrieving the rankings for a particular set of genes.
RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings WHERE geneID IN ({0:s}) ORDER BY geneID;"
# SQL query that retrieves the ordered list of features in the database.
FEATURE_IDS_QUERY = r"SELECT motifName FROM motifs ORDER BY idx;"


class RankingDatabase:
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
        assert os.path.exists(fname), "Database {0:s} doesn't exist.".format(fname)
        assert name, "Name must be specified."
        assert nomenclature, "Nomenclature must be specified."

        self._fname = fname
        self._name = name
        self._nomenclature = nomenclature

        # Read-only view on SQLite database.
        self._uri = 'file:{}?mode=ro'.format(os.path.abspath(fname))

        def fetch_features():
            with sqlite3.connect(self._uri, uri=True) as db:
                cursor = db.cursor()
                features = tuple(map(itemgetter(0), cursor.execute(FEATURE_IDS_QUERY).fetchall()))
                cursor.close()
            return features
        self._features = fetch_features()

        def fetch_gene_count():
            with sqlite3.connect(self._uri, uri=True) as db:
                cursor = db.cursor()
                count = cursor.execute(GENE_ID_COUNT_QUERY).fetchone()
                cursor.close()
            return count
        self._gene_count = fetch_gene_count()

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
    def features(self) -> Tuple[str]:
        """
        List of regulatory features for which whole genome rankings are available in this database.
        """
        return self._features

    def load(self, gs: GeneSignature) -> (np.ndarray,np.ndarray,np.ndarray):
        """
        Load the ranking of the genes in the supplied signature for all features in this database.

        :param gs: The gene signature.
        :return: A tuple of numpy ndarrays: the vector of features, the vector of genes and the (n_features, n_genes)
            matrix containing the actual rankings of these genes for each regulatory feature in the database.
        """
        assert gs, "A gene signature must be supplied"

        def quoted_csv(values):
            # Escape single quotes (') by using (''), because sometimes ID's contain a single quote.
            def quote(value):
                return "'" + value.replace("'", "''") + "'"
            return ','.join(map(quote, values))

        # Pre-allocate the matrix.
        rankings = np.empty(shape=(len(self._features), len(gs)), dtype=self._dtype)
        with sqlite3.connect(self._uri, uri=True) as db:
            cursor = db.cursor()
            cursor.execute(RANKINGS_QUERY.format(quoted_csv(gs.genes)))
            row_idx = 0
            genes = []
            for gene, ranking in cursor:
                rankings[row_idx, :] = np.frombuffer(ranking, dtype=self._dtype)
                genes.append(gene)
                row_idx += 1
            cursor.close()

        return np.array(self._features, dtype='|S255'), np.array(genes, dtype='|S255'), rankings

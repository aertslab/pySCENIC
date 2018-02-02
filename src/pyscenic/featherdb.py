# -*- coding: utf-8 -*-

import os
import pandas as pd
from feather.api import write_dataframe, FeatherReader
from typing import Type, Tuple
from .genesig import GeneSignature
from .sqlitedb import SQLiteRankingDatabase
from .rnkdb import RankingDatabase


def convert2feather(fname: str, out_folder: str, name: str, nomenclature: str, extension: str="feather") -> str:
    """
    Convert a whole genome rankings database to a feather format based database.

    More information on this format can be found here:
    .. feather-format: https://blog.rstudio.com/2016/03/29/feather/

    :param fname: The filename of the legacy
    :param out_folder: The name of the folder to write the new database to.
    :param name: The name of the rankings database.
    :param nomenclature: The nomenclature used for the genes.
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
    db = SQLiteRankingDatabase(fname=fname, name=name, nomenclature=nomenclature)
    features, genes, rankings = db.load_full()
    # Genes must be columns because feather is a column-oriented format.
    df = pd.DataFrame(index=features, columns=genes, data=rankings)
    # Nomenclature is stored as the name of the column-index of the dataframe.
    df.columns.name = nomenclature
    # The name of the database of rankings is stored as the name of the index of the dataframe.
    df.index.name = name
    write_dataframe(df, feather_fname)
    return feather_fname


class FeatherRankingDatabase(RankingDatabase):
    def __init__(self, fname: str, name: str = None, nomenclature: str = None):
        """
        Create a new feather database.

        :param fname: The filename of the database.
        :param name: The name of the database.
        :param nomenclature: The nomenclature used for the genes that are being ranked.
        """
        super().__init__(fname, name=name, nomenclature=nomenclature)
        self._reader = FeatherReader(fname)

    @property
    def total_genes(self) -> int:
        return self._reader.num_columns

    @property
    def genes(self) -> Tuple[str]:
        # noinspection PyTypeChecker
        return tuple(self._reader.get_column_name(idx) for idx in range(self.total_genes))

    def load_full(self) -> pd.DataFrame:
        return self._reader.read()

    def load(self, gs: Type[GeneSignature]) -> pd.DataFrame:
        return self._reader.read(columns=gs.genes)

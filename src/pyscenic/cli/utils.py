# -*- coding: utf-8 -*-

import os
import pickle
import pandas as pd
import loompy as lp
from typing import Type, List, Sequence, Tuple
from pyscenic.genesig import GeneSignature
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs, load_from_yaml

__all__ = ['save_matrix', 'load_exp_matrix', 'load_signatures']


ATTRIBUTE_NAME_CELL_IDENTIFIER = "CellID"
ATTRIBUTE_NAME_GENE = "Gene"


def save_df_as_loom(df: pd.DataFrame, fname: str) -> None:
    """
    Save pandas dataframe as single layer loom file. Can be used to save expression matrix or AUC value matrix
    as binary loom file.

    :param df: The 2-dimensional dataframe (rows = cells x columns = genes).
    :param fname: The name of the loom file to create.
    """
    assert df.ndim == 2
    # The orientation of the loom file is always:
    #   - Columns represent cells or aggregates of cells
    # 	- Rows represent genes
    column_attrs = { ATTRIBUTE_NAME_CELL_IDENTIFIER: df.index.values.astype('str'), }
    row_attrs = { ATTRIBUTE_NAME_GENE: df.columns.values.astype('str'), }
    lp.create(filename=fname,
              layers=df.T.values,
              row_attrs=row_attrs,
              col_attrs=column_attrs)


def load_exp_matrix_as_loom(fname) -> pd.DataFrame:
    """
    Load expression matrix from loom file.

    :param fname: The name of the loom file to load.
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    with lp.connect(fname) as ds:
        # The orientation of the loom file is always:
        #   - Columns represent cells or aggregates of cells
        # 	- Rows represent genes
        return pd.DataFrame(data=ds[:, :],
                            index=ds.ra[ATTRIBUTE_NAME_GENE],
                            columns=ds.ca[ATTRIBUTE_NAME_CELL_IDENTIFIER]).T


FILE_EXTENSION2SEPARATOR = {
    '.tsv': '\t',
    '.csv': ','
}


def load_exp_matrix(fname: str, transpose: bool = False) -> pd.DataFrame:
    """
    Load expression matrix from disk.

    Supported file formats are CSV, TSV and LOOM.

    :param fname: The name of the file that contains the expression matrix.
    :param transpose: Is the expression matrix stored as (rows = genes x columns = cells)?
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    extension = os.path.splitext(fname)[1].lower()
    if extension in FILE_EXTENSION2SEPARATOR.keys():
        df = pd.read_csv(fname, sep=FILE_EXTENSION2SEPARATOR[extension], header=0, index_col=0)
        return df.T if transpose else df
    elif extension == '.loom':
        return load_exp_matrix_as_loom(fname)
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))


def save_matrix(df: pd.DataFrame, fname: str, transpose: bool = False) -> None:
    """
    Save matrix to disk.

    Supported file formats are CSV, TSV and LOOM.

    :param df: A 2-dimensional dataframe (rows = cells x columns = genes).
    :param fname: The name of the file to be written.
    :param transpose: Should the expression matrix be stored as (rows = genes x columns = cells)?
    """
    extension = os.path.splitext(fname)[1].lower()
    if extension in FILE_EXTENSION2SEPARATOR.keys():
        (df.T if transpose else df).to_csv(fname, sep=FILE_EXTENSION2SEPARATOR[extension])
    elif extension == '.loom':
        return save_df_as_loom(df, fname)
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))


LINE_LIMIT = 100


def guess_separator(fname: str) -> str:
    with open(fname, 'r') as f:
        lines = f.readline(LINE_LIMIT)

    def count_columns(sep):
        return [len(line.split(sep)) for line in lines if not line.strip().startswith('#') and line.strip()]

    # Check if '\t' is used:
    for sep in ('\t', ';', ','):
        if min(count_columns(sep)) >= 3:
            return sep
    return ValueError("Unknown file format \"{}\".".format(fname))


def load_signatures(fname: str) -> Sequence[Type[GeneSignature]]:
    """
    Load genes signatures from disk.

    Supported file formats are GMT, DAT (pickled), YAML or CSV (enriched motifs).

    :param fname: The name of the file that contains the signatures.
    :return: A list of gene signatures.
    """
    extension = os.path.splitext(fname)[1].lower()
    if extension in FILE_EXTENSION2SEPARATOR.keys():
        return df2regulons(load_motifs(fname, sep=FILE_EXTENSION2SEPARATOR[extension]))
    elif extension in {'.yaml', '.yml'}:
        return load_from_yaml(fname)
    elif extension == '.gmt':
        sep = guess_separator(fname)
        return GeneSignature.from_gmt(fname,
                                  field_separator=sep,
                                  gene_separator=sep)
    elif extension == '.dat':
        with open(fname, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))
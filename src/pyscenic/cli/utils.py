# -*- coding: utf-8 -*-

import pandas as pd
import loompy as lp


ATTRIBUTE_NAME_CELL_IDENTIFIER = "CellID"
ATTRIBUTE_NAME_GENE = "Gene"


def save_df_as_loom(df: pd.DataFrame, fname: str) -> None:
    """
    Save pandas dataframe as single layer loom file. Can be used to save expression matrix or AUC value matrix
    as binary loom file.

    :param df: The 2-dimensional dataframe (rows = cells x columns = genes).
    :param fname: The name of the loom file to create.
    """
    assert pd.ndim == 2
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

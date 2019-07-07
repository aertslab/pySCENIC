# -*- coding: utf-8 -*-

import os
import pickle
import json
import zlib
import base64
import numpy as np
import pandas as pd
import loompy as lp
from operator import attrgetter
from typing import Type, Sequence
from pyscenic.genesig import GeneSignature
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs, load_from_yaml, save_to_yaml
from pyscenic.binarization import binarize


__all__ = ['save_matrix', 'load_exp_matrix', 'load_signatures', 'save_enriched_motifs', 'load_adjacencies',
           'load_modules', 'append_auc_mtx']


ATTRIBUTE_NAME_CELL_IDENTIFIER = "CellID"
ATTRIBUTE_NAME_GENE = "Gene"
ATTRIBUTE_NAME_REGULONS_AUC = "RegulonsAUC"
ATTRIBUTE_NAME_REGULONS = "Regulons"
ATTRIBUTE_NAME_METADATA = "MetaData"


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


def load_exp_matrix_as_loom(fname,
                            attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
                            attribute_name_gene: str = ATTRIBUTE_NAME_GENE) -> pd.DataFrame:
    """
    Load expression matrix from loom file.

    :param fname: The name of the loom file to load.
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    with lp.connect(fname,mode='r',validate=False) as ds:
        # The orientation of the loom file is always:
        #   - Columns represent cells or aggregates of cells
        # 	- Rows represent genes
        return pd.DataFrame(data=ds[:, :],
                            index=ds.ra[attribute_name_gene],
                            columns=ds.ca[attribute_name_cell_id]).T


FILE_EXTENSION2SEPARATOR = {
    '.tsv': '\t',
    '.csv': ','
}


def load_exp_matrix(fname: str, transpose: bool = False,
                    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
                    attribute_name_gene: str = ATTRIBUTE_NAME_GENE) -> pd.DataFrame:
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
        return load_exp_matrix_as_loom(fname, attribute_name_cell_id, attribute_name_gene)
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


def guess_separator(fname: str) -> str:
    with open(fname, 'r') as f:
        lines = f.readlines()

    def count_columns(sep):
        return [len(line.split(sep)) for line in lines if not line.strip().startswith('#') and line.strip()]

    # Check if '\t' is used:
    for sep in ('\t', ';', ','):
        if min(count_columns(sep)) >= 3:
            return sep
    raise ValueError("Unknown file format \"{}\".".format(fname))


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
    elif extension.endswith('.gmt'):
        sep = guess_separator(fname)
        return GeneSignature.from_gmt(fname,
                                  field_separator=sep,
                                  gene_separator=sep)
    elif extension == '.dat':
        with open(fname, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))


def save_enriched_motifs(df, fname:str) -> None:
    """
    Save enriched motifs.

    Supported file formats are CSV, TSV, GMT, DAT (pickle), JSON or YAML.

    :param df:
    :param fname:
    :return:
    """
    extension = os.path.splitext(fname)[1].lower()
    if extension in FILE_EXTENSION2SEPARATOR.keys():
        df.to_csv(fname, sep=FILE_EXTENSION2SEPARATOR[extension])
    else:
        regulons = df2regulons(df)
        if extension == '.json':
            name2targets = {r.name: list(r.gene2weight.keys()) for r in regulons}
            with open(fname, 'w') as f:
                f.write(json.dumps(name2targets))
        elif extension == '.dat':
            pickle.dump(regulons, fname)
        elif extension == '.gmt':
            GeneSignature.to_gmt(fname, regulons)
        elif extension in {'.yaml', '.yml'}:
            save_to_yaml(regulons, fname)
        else:
            raise ValueError("Unknown file format \"{}\".".format(fname))


def load_adjacencies(fname: str) -> pd.DataFrame:
    extension = os.path.splitext(fname)[1].lower().lower()
    return pd.read_csv(fname, sep=FILE_EXTENSION2SEPARATOR[extension], dtype={0:str,1:str,2:np.float64}, keep_default_na=False )


def load_modules(fname: str) -> Sequence[Type[GeneSignature]]:
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    if fname.endswith('.yaml') or fname.endswith('.yml'):
        return load_from_yaml(fname)
    elif fname.endswith('.dat'):
        with open(fname, 'rb') as f:
            return pickle.load(f)
    elif fname.endswith('.gmt'):
        sep = guess_separator(fname)
        return GeneSignature.from_gmt(fname,
                                      field_separator=sep,
                                      gene_separator=sep)
    else:
        raise ValueError("Unknown file format for \"{}\".".format(fname))


def decompress_meta(meta):
    try:
        meta = meta.decode('ascii')
        return json.loads(zlib.decompress(base64.b64decode(meta)))
    except AttributeError:
        return json.loads(zlib.decompress(base64.b64decode(meta.encode('ascii'))).decode('ascii'))


def compress_meta(meta):
    return base64.b64encode(zlib.compress(json.dumps(meta).encode('ascii'))).decode('ascii')


def append_auc_mtx(fname: str, auc_mtx: pd.DataFrame, regulons: Sequence[Type[GeneSignature]]) -> None:
    """

    Append AUC matrix to loom file.

    :param fname: The name of loom file to be append to.
    :param auc_mtx: The matrix that contains the AUC values.
    :param regulons: Collection of regulons that were used for calculation of the AUC values.
    """
    # Fetch sequence logo from regulon's context.
    def fetch_logo(context):
        for elem in context:
            if elem.endswith('.png'):
                return elem
        return ""
    try:
        name2logo = {reg.name: fetch_logo(reg.context) for reg in regulons}
    except AttributeError:
        name2logo = {}

    # Binarize matrix for AUC thresholds.
    _, auc_thresholds = binarize(auc_mtx)
    regulon_thresholds = [{"regulon": name,
                           "defaultThresholdValue":(threshold if isinstance(threshold, float) else threshold[0]),
                           "defaultThresholdName": "guassian_mixture_split",
                           "allThresholds": {"guassian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])},
                           "motifData": name2logo.get(name, "")} for name, threshold in auc_thresholds.iteritems()]

    # Calculate the number of genes per cell.
    ex_mtx = load_exp_matrix(fname)
    binary_mtx = ex_mtx.copy()
    binary_mtx[binary_mtx != 0] = 1.0
    ngenes = binary_mtx.sum(axis=1).astype(int)

    # Encode genes in regulons as "binary" membership matrix.
    genes = np.array(ex_mtx.columns)
    n_genes = len(genes)
    n_regulons = len(regulons)
    data = np.zeros(shape=(n_genes, n_regulons), dtype=int)
    for idx, regulon in enumerate(regulons):
        data[:, idx] = np.isin(genes, regulon.genes).astype(int)
    regulon_assignment = pd.DataFrame(data=data,
                                      index=ex_mtx.columns,
                                      columns=list(map(attrgetter('name'), regulons)))

    # Create meta-data structure.
    def create_structure_array(df):
        # Create a numpy structured array
        return np.array([tuple(row) for row in df.values],
                        dtype=np.dtype(list(zip(df.columns, df.dtypes))))

    with lp.connect(fname, validate=False) as ds:
        # The orientation of the loom file is always:
        #   - Columns represent cells or aggregates of cells
        # 	- Rows represent genes
        ds.ca[ATTRIBUTE_NAME_REGULONS_AUC] = create_structure_array(auc_mtx)
        ds.ra[ATTRIBUTE_NAME_REGULONS] = create_structure_array(regulon_assignment)
        if ATTRIBUTE_NAME_METADATA in ds.attrs:
            try:
                meta_data = json.loads(ds.attrs[ATTRIBUTE_NAME_METADATA])
            except json.decoder.JSONDecodeError:
                meta_data = decompress_meta(ds.attrs[ATTRIBUTE_NAME_METADATA])
        else:
            meta_data = {}
        meta_data["regulonThresholds"] = regulon_thresholds
        ds.attrs[ATTRIBUTE_NAME_METADATA] = compress_meta(meta_data)

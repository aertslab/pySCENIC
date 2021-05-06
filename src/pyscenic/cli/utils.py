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
from ctxcore.genesig import GeneSignature, openfile
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs, load_from_yaml, save_to_yaml
from pyscenic.binarization import binarize
from pathlib import Path, PurePath


__all__ = [
    'save_matrix',
    'load_exp_matrix',
    'load_signatures',
    'save_enriched_motifs',
    'load_adjacencies',
    'load_modules',
    'append_auc_mtx',
]


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
    column_attrs = {
        ATTRIBUTE_NAME_CELL_IDENTIFIER: df.index.values.astype('str'),
    }
    row_attrs = {
        ATTRIBUTE_NAME_GENE: df.columns.values.astype('str'),
    }
    lp.create(filename=fname, layers=df.T.values, row_attrs=row_attrs, col_attrs=column_attrs)


def load_exp_matrix_as_loom(
    fname,
    return_sparse=False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from loom file.

    :param fname: The name of the loom file to load.
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    if return_sparse:
        with lp.connect(fname, mode='r', validate=False) as ds:
            ex_mtx = ds.layers[''].sparse().T.tocsc()
            genes = pd.Series(ds.ra[attribute_name_gene])
            cells = ds.ca[attribute_name_cell_id]
        return ex_mtx, genes, cells

    else:
        with lp.connect(fname, mode='r', validate=False) as ds:
            # The orientation of the loom file is always:
            #   - Columns represent cells or aggregates of cells
            # 	- Rows represent genes
            return pd.DataFrame(
                data=ds[:, :], index=ds.ra[attribute_name_gene], columns=ds.ca[attribute_name_cell_id]
            ).T


def suffixes_to_separator(extension):
    if '.csv' in extension:
        return ','
    if '.tsv' in extension:
        return '\t'


def is_valid_suffix(extension, method):
    assert isinstance(extension, list), 'extension should be of type "list"'
    if method in ['grn', 'aucell']:
        valid_extensions = ['.csv', '.tsv', '.loom', '.h5ad']
    elif method == 'ctx':
        valid_extensions = ['.csv', '.tsv']
    elif method == 'ctx_yaml':
        valid_extensions = ['.yaml', '.yml']
    if len(set(extension).intersection(valid_extensions)) > 0:
        return True
    else:
        return False


def load_exp_matrix(
    fname: str,
    transpose: bool = False,
    return_sparse: bool = False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from disk.

    Supported file formats are CSV, TSV and LOOM.

    :param fname: The name of the file that contains the expression matrix.
    :param transpose: Is the expression matrix stored as (rows = genes x columns = cells)?
    :param return_sparse: Returns a sparse matrix when loading from loom
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    extension = PurePath(fname).suffixes
    if is_valid_suffix(extension, 'grn'):
        if '.loom' in extension:
            return load_exp_matrix_as_loom(fname, return_sparse, attribute_name_cell_id, attribute_name_gene)
        elif '.h5ad' in extension:
            from anndata import read_h5ad

            adata = read_h5ad(filename=fname, backed='r')
            if return_sparse:
                # expr, gene, cell:
                return adata.X.value, adata.var_names.values, adata.obs_names.values
            else:
                return pd.DataFrame(
                    adata.X.value.todense(), index=adata.obs_names.values, columns=adata.var_names.values
                )

        else:
            df = pd.read_csv(fname, sep=suffixes_to_separator(extension), header=0, index_col=0)
            return df.T if transpose else df
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
    extension = PurePath(fname).suffixes
    if is_valid_suffix(extension, 'aucell'):
        if '.loom' in extension:
            return save_df_as_loom(df, fname)
        else:
            (df.T if transpose else df).to_csv(fname, sep=suffixes_to_separator(extension))
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))


def guess_separator(fname: str) -> str:
    with openfile(fname, 'r') as f:
        lines = f.readlines()

    # decode if gzipped file:
    for i, x in enumerate(lines):
        if isinstance(x, (bytes, bytearray)):
            lines[i] = x.decode()

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
    extension = PurePath(fname).suffixes
    if is_valid_suffix(extension, 'ctx'):
        # csv/tsv
        return df2regulons(load_motifs(fname, sep=suffixes_to_separator(extension)))
    elif is_valid_suffix(extension, 'ctx_yaml'):
        return load_from_yaml(fname)
    elif '.gmt' in extension:
        sep = guess_separator(fname)
        return GeneSignature.from_gmt(fname, field_separator=sep, gene_separator=sep)
    elif '.dat' in extension:
        with openfile(fname, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))


def save_enriched_motifs(df, fname: str) -> None:
    """
    Save enriched motifs.

    Supported file formats are CSV, TSV, GMT, DAT (pickle), JSON or YAML.

    :param df:
    :param fname:
    :return:
    """
    extension = PurePath(fname).suffixes
    if is_valid_suffix(extension, 'ctx'):
        df.to_csv(fname, sep=suffixes_to_separator(extension))
    else:
        regulons = df2regulons(df)
        if '.json' in extension:
            name2targets = {r.name: list(r.gene2weight.keys()) for r in regulons}
            with openfile(fname, 'w') as f:
                f.write(json.dumps(name2targets))
        elif '.dat' in extension:
            with openfile(fname, 'wb') as f:
                pickle.dump(regulons, f)
        elif '.gmt' in extension:
            GeneSignature.to_gmt(fname, regulons)
        elif is_valid_suffix(extension, 'ctx_yaml'):
            save_to_yaml(regulons, fname)
        else:
            raise ValueError("Unknown file format \"{}\".".format(fname))


def load_adjacencies(fname: str) -> pd.DataFrame:
    extension = PurePath(fname).suffixes
    return pd.read_csv(
        fname, sep=suffixes_to_separator(extension), dtype={0: str, 1: str, 2: np.float64}, keep_default_na=False
    )


def load_modules(fname: str) -> Sequence[Type[GeneSignature]]:
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    extension = PurePath(fname).suffixes
    if is_valid_suffix(extension, 'ctx_yaml'):
        return load_from_yaml(fname)
    elif '.dat' in extension:
        with openfile(fname, 'rb') as f:
            return pickle.load(f)
    elif '.gmt' in extension:
        return GeneSignature.from_gmt(fname)
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


def append_auc_mtx(
    fname: str,
    ex_mtx: pd.DataFrame,
    auc_mtx: pd.DataFrame,
    regulons: Sequence[Type[GeneSignature]],
    seed=None,
    num_workers=1,
) -> None:
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
    _, auc_thresholds = binarize(auc_mtx, seed=seed, num_workers=num_workers)
    regulon_thresholds = [
        {
            "regulon": name,
            "defaultThresholdValue": (threshold if isinstance(threshold, float) else threshold[0]),
            "defaultThresholdName": "gaussian_mixture_split",
            "allThresholds": {"gaussian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])},
            "motifData": name2logo.get(name, ""),
        }
        for name, threshold in auc_thresholds.iteritems()
    ]

    # Calculate the number of genes per cell.
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
    regulon_assignment = pd.DataFrame(data=data, index=ex_mtx.columns, columns=list(map(attrgetter('name'), regulons)))

    # Create meta-data structure.
    def create_structure_array(df):
        # Create a numpy structured array
        return np.array([tuple(row) for row in df.values], dtype=np.dtype(list(zip(df.columns, df.dtypes))))

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

# coding=utf-8

import os
import numpy as np
import pandas as pd
import loompy as lp
from umap import UMAP
from .aucell import aucell
from .genesig import Regulon
from typing import List, Mapping
from operator import attrgetter
from multiprocessing import cpu_count


def export2loom(ex_mtx: pd.DataFrame, regulons: List[Regulon], cell_annotations: Mapping[str,str],
                out_fname: str, num_cores=cpu_count()):
    """
    Create a loom file for a single cell experiment to be used in SCope.

    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param regulons: A list of Regulons.
    :param cell_annotations: A dictionary that maps a cell ID to its corresponding cell type annotation.
    :param out_fname: The name of the file to create.
    :param num_cores: The number of cores to use for AUCell regulon enrichment.
    """
    # Information on the general loom file format: http://linnarssonlab.org/loompy/format/index.html
    # Information on the SCope specific alterations: https://github.com/aertslab/SCope/wiki/Data-Format

    # TODO: Add regulon thresholds: not mandatory!

    # Calculate regulon enrichment per cell using AUCell.
    auc_mtx = aucell(ex_mtx, regulons, num_cores) # (n_cells x n_regulons)

    # Create an embedding based on UMAP (similar to tSNE but faster).
    umap_embedding_mtx = UMAP().fit_transform(auc_mtx) # (n_cells, 2)

    # Calculate the number of genes per cell.
    binary_mtx = ex_mtx.copy()
    binary_mtx[binary_mtx != 0] = 1.0
    ngenes = binary_mtx.sum(axis=1)

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
        return np.array([tuple(row) for row in df.as_matrix()],
                        dtype=np.dtype(list(zip(df.columns, df.dtypes))))


    nomenclatures = set(map(attrgetter('nomenclature'), regulons))
    assert len(nomenclatures) == 1

    title = os.path.splitext(os.path.basename(out_fname))[0]

    column_attrs = {
        "CellID": list(ex_mtx.index),
        "nGene": list(ngenes.values),
        "Embedding": create_structure_array(umap_embedding_mtx.T),
        "RegulonsAUC": create_structure_array(auc_mtx),
        "ClusterID": list(ex_mtx.index.rename(cell_annotations))
        }
    row_attrs = {
        "Gene": list(ex_mtx.columns),
        "Regulons": regulon_assignment,
        }
    general_attrs = {
        "title": title,
        "MetaData": {
            "embeddings": [{
                "id": 0,
                "name": "UMAP (default)",
            }]},
        "Genome": next(iter(nomenclatures))}

    # Create loom file for use with the SCope tool.
    # The loom file format opted for rows as genes to facilitate growth along the column axis (i.e add more cells)
    # PySCENIC chose a different orientation because of limitation set by the feather format: selectively reading
    # information from disk can only be achieved via column selection. For the ranking databases this is of utmost
    # importance.
    fh = lp.create(filename=out_fname,
              layers=ex_mtx.T.values,
              row_attrs=row_attrs,
              col_attrs=column_attrs,
              file_attrs=general_attrs)
    fh.close()

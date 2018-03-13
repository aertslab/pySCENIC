# -*- coding: utf-8 -*-

import pandas as pd
from .recovery import enrichment4cells
from tqdm import tqdm
from typing import Sequence, Type
from .genesig import GeneSignature
from joblib import Parallel, delayed, dump, load
from multiprocessing import cpu_count
import tempfile
import os
import logging


LOGGER = logging.getLogger(__name__)


def create_rankings(ex_mtx: pd.DataFrame) -> pd.DataFrame:
    """
    Create a whole genome rankings dataframe from a single cell expression profile dataframe.

    :param ex_mtx: The expression profile matrix. The rows should correspond to different cells, the columns to different
        genes.
    :return: A genome rankings dataframe.
    """
    return ex_mtx.rank(axis=1, ascending=False, method='first').astype('int64')


enrichment = enrichment4cells


def _enrichment(rnk_mtx_as_ndarray, module, genes, features, rank_threshold, auc_threshold):
    return enrichment4cells(pd.DataFrame(data=rnk_mtx_as_ndarray, columns=genes, index=features),
                     module, rank_threshold, auc_threshold)


def aucell(exp_mtx: pd.DataFrame, modules: Sequence[Type[GeneSignature]],
           rank_threshold: int = 5000, auc_threshold: float = 0.05,
           noweights: bool = False, num_cores: int = cpu_count()) -> pd.DataFrame:
    """
    Calculate enrichment of regulomes for single cells.

    :param exp_mtx: The expression matrix (n_cells x n_genes).
    :param modules: The regulomes.
    :param rank_threshold:
    :param auc_threshold:
    :param noweights: Should the weights of the genes part of regulome be used in calculation of enrichment?
    :param num_cores:
    :return:
    """
    df_rnk = create_rankings(exp_mtx)

    if num_cores == 1:
        # Show progress bar ...
        aucs = pd.concat([enrichment(df_rnk, module if noweights else module.noweights(),
                                 rank_threshold=rank_threshold,
                                 auc_threshold=auc_threshold) for module in tqdm(modules)])
    else:
        temp_fname = tempfile.mktemp()
        try:
            # Decompose rankings dataframe so the core rankings ndarray can be shared amongst
            # processes using a memory mapped file.
            genes = df_rnk.columns.values
            features = df_rnk.index.values
            rnk_mtx = df_rnk.as_matrix()

            # Dump the input data to disk to free the memory
            dump(rnk_mtx, temp_fname)

            # Release the reference on the original in memory array and replace it
            # by a reference to the memmap array so that the garbage collector can
            # release the memory before forking. gc.collect() is internally called
            # in Parallel just before forking.
            rnk_mtx = load(temp_fname, mmap_mode='r')

            if noweights:
                modules = list(map(lambda m: m.noweights(), modules))
            aucs = pd.concat(Parallel(n_jobs=num_cores)(delayed(_enrichment)
                                      (rnk_mtx, module,
                                            genes, features, rank_threshold, auc_threshold) for module in tqdm(modules)))
        finally:
            try:
                os.remove(temp_fname)
            except:
                logging.warning("Failed to delete \"{}\".".format(temp_fname))

    return aucs.unstack("Regulome")

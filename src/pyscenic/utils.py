# -*- coding: utf-8 -*-

import pandas as pd
from .genesig import Regulome, GeneSignature
from collections import defaultdict, Counter
from itertools import chain
import numpy as np
from functools import partial
from typing import Sequence, Type
import yaml


COLUMN_NAME_MOTIF_SIMILARITY_QVALUE = 'motif_similarity_qvalue'
COLUMN_NAME_ORTHOLOGOUS_IDENTITY = 'orthologous_identity'


def load_motif_annotations(fname: str,
                           column_names=('#motif_id', 'gene_name',
                                         COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
                                         'description'),
                           motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0) -> pd.DataFrame:
    """
    Load motif annotations from a motif2TF snapshot.

    :param fname: the snapshot taken from motif2TF.
    :param column_names: the names of the columns in the snapshot to load.
    :param motif_similarity_fdr: The maximum False Discovery Rate to find factor annotations for enriched motifs.
    :param orthologuous_identity_threshold: The minimum orthologuous identity to find factor annotations
        for enriched motifs.
    :return: A dataframe.
    """
    # Create a MultiIndex for the index combining unique gene name and motif ID. This should facilitate
    # later merging.
    df = pd.read_csv(fname, sep='\t', index_col=[1,0], usecols=column_names)
    df = df[(df[COLUMN_NAME_MOTIF_SIMILARITY_QVALUE] <= motif_similarity_fdr) &
            (df[COLUMN_NAME_ORTHOLOGOUS_IDENTITY] >= orthologuous_identity_threshold)]
    return df


COLUMN_NAME_TF = "TF"
COLUMN_NAME_TARGET = "target"
COLUMN_NAME_WEIGHT = "importance"
COLUMN_NAME_CORRELATION = "correlation"
RHO_THRESHOLD = 0.3


def add_correlation(adjacencies: pd.DataFrame, ex_mtx: pd.DataFrame) -> pd.DataFrame:
    """
    Add correlation in expression levels between target and factor.

    :param adjacencies: The dataframe with the TF-target links.
    :param ex_mtx: The expression matrix (n_genes x n_cells).
    :return: The adjacencies dataframe with an extra column.
    """

    # Calculate Pearson correlation to infer repression or activation.
    corr_mtx = pd.DataFrame(index=ex_mtx.index, columns=ex_mtx.index, data=np.corrcoef(ex_mtx.values))

    # Add "correlation" column to adjacencies dataframe.
    def add_regulation(row, corr_mtx):
        tf = row[COLUMN_NAME_TF]
        target = row[COLUMN_NAME_TARGET]
        rho = corr_mtx[tf][target]
        return int(rho > RHO_THRESHOLD) - int(rho < RHO_THRESHOLD)

    adjacencies[COLUMN_NAME_CORRELATION] = adjacencies.apply(partial(add_regulation, corr_mtx=corr_mtx), axis=1)

    return adjacencies


def modules4thr(adjacencies, threshold, nomenclature="MGI"):
    """

    :param adjacencies:
    :param threshold:
    :param nomenclature:
    :return:
    """
    for tf_name, df_grp in adjacencies[adjacencies[COLUMN_NAME_WEIGHT] > threshold].groupby(by=COLUMN_NAME_TF):
        if len(df_grp) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=frozenset(["weight>{}".format(threshold)]),
                transcription_factor=tf_name,
                gene2weights=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


def modules4top_targets(adjacencies, n, nomenclature="MGI"):
    """

    :param adjacencies:
    :param n:
    :param nomenclature:
    :return:
    """
    for tf_name, df_grp in adjacencies.groupby(by=COLUMN_NAME_TF):
        module = df_grp.nlargest(n, COLUMN_NAME_WEIGHT)
        if len(module) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=frozenset(["top{}".format(n)]),
                transcription_factor=tf_name,
                gene2weights=list(zip(module[COLUMN_NAME_TARGET].values, module[COLUMN_NAME_WEIGHT].values)))


def modules4top_factors(adjacencies, n, nomenclature="MGI"):
    """

    :param adjacencies:
    :param n:
    :param nomenclature:
    :return:
    """
    df = adjacencies.groupby(by=COLUMN_NAME_TARGET).apply(lambda grp: grp.nlargest(n, COLUMN_NAME_WEIGHT))
    for tf_name, df_grp in df.groupby(by=COLUMN_NAME_TF):
        if len(df_grp) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=frozenset(["top{}perTarget".format(n)]),
                transcription_factor=tf_name,
                gene2weights=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


def modules_from_genie3(adjacencies: pd.DataFrame, nomenclature: str,
                        thresholds=(0.001,0.005),
                        top_n_targets=(50,),
                        top_n_regulators=(5,10,50)):
    """
    
    :param adjacencies:
    :param nomenclature:
    :param thresholds:
    :param top_n_targets:
    :param top_n_regulators:
    :return:
    """
    yield from chain(chain.from_iterable(modules4thr(adjacencies, thr, nomenclature) for thr in thresholds),
                     chain.from_iterable(modules4top_targets(adjacencies, n, nomenclature) for n in top_n_targets),
                     chain.from_iterable(modules4top_factors(adjacencies, n, nomenclature) for n in top_n_regulators))


def save_as_yaml(signatures: Sequence[Type[GeneSignature]], fname: str):
    """

    :param signatures:
    :return:
    """
    with open(fname, 'w') as f:
        f.write(yaml.dump(signatures))


def load_from_yaml(fname: str) -> Sequence[Type[GeneSignature]]:
    """

    :param fname:
    :return:
    """
    with open(fname, 'r') as f:
        return yaml.load(f.read())

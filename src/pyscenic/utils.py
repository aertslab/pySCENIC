# -*- coding: utf-8 -*-

import pandas as pd
from urllib.parse import urljoin
from .genesig import Regulome, GeneSignature
from itertools import chain
import numpy as np
from functools import partial
from typing import Sequence, Type
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


COLUMN_NAME_TF = "TF"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_MOTIF_URL = "MotifURL"
COLUMN_NAME_MOTIF_SIMILARITY_QVALUE = 'MotifSimilarityQvalue'
COLUMN_NAME_ORTHOLOGOUS_IDENTITY = 'OrthologousIdentity'
COLUMN_NAME_ANNOTATION = 'Annotation'


def load_motif_annotations(fname: str,
                           column_names=('#motif_id', 'gene_name',
                                         'motif_similarity_qvalue', 'orthologous_identity', 'description'),
                           motif_similarity_fdr: float = 0.001,
                           orthologous_identity_threshold: float = 0.0) -> pd.DataFrame:
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
    df.index.names = [COLUMN_NAME_TF, COLUMN_NAME_MOTIF_ID]
    df.rename(columns={'motif_similarity_qvalue': COLUMN_NAME_MOTIF_SIMILARITY_QVALUE,
                       'orthologous_identity': COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
                       'description': COLUMN_NAME_ANNOTATION }, inplace=True)
    df = df[(df[COLUMN_NAME_MOTIF_SIMILARITY_QVALUE] <= motif_similarity_fdr) &
            (df[COLUMN_NAME_ORTHOLOGOUS_IDENTITY] >= orthologous_identity_threshold)]
    return df


COLUMN_NAME_TARGET = "target"
COLUMN_NAME_WEIGHT = "importance"
COLUMN_NAME_CORRELATION = "correlation"
RHO_THRESHOLD = 0.03


def add_correlation(adjacencies: pd.DataFrame, ex_mtx: pd.DataFrame, rho_threshold=RHO_THRESHOLD) -> pd.DataFrame:
    """
    Add correlation in expression levels between target and factor.

    :param adjacencies: The dataframe with the TF-target links.
    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param rho_threshold: The threshold on the correlation to decide if a target gene is activated
        (rho > `rho_threshold`) or repressed (rho < -`rho_threshold`).
    :return: The adjacencies dataframe with an extra column.
    """
    assert rho_threshold > 0, "rho_threshold should be greater than 0."

    # Calculate Pearson correlation to infer repression or activation.
    corr_mtx = pd.DataFrame(index=ex_mtx.columns, columns=ex_mtx.columns, data=np.corrcoef(ex_mtx.values.T))

    # Add "correlation" column to adjacencies dataframe.
    def add_regulation(row, corr_mtx):
        tf = row[COLUMN_NAME_TF]
        target = row[COLUMN_NAME_TARGET]
        rho = corr_mtx[tf][target]
        return int(rho > rho_threshold) - int(rho < -rho_threshold)

    adjacencies[COLUMN_NAME_CORRELATION] = adjacencies.apply(partial(add_regulation, corr_mtx=corr_mtx), axis=1)

    return adjacencies


def modules4thr(adjacencies, threshold, nomenclature="MGI", context=frozenset()):
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
                context=frozenset(["weight>{}".format(threshold)]).union(context),
                transcription_factor=tf_name,
                gene2weights=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


def modules4top_targets(adjacencies, n, nomenclature="MGI", context=frozenset()):
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
                context=frozenset(["top{}".format(n)]).union(context),
                transcription_factor=tf_name,
                gene2weights=list(zip(module[COLUMN_NAME_TARGET].values, module[COLUMN_NAME_WEIGHT].values)))


def modules4top_factors(adjacencies, n, nomenclature="MGI", context=frozenset()):
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
                name=tf_name,
                nomenclature=nomenclature,
                context=frozenset(["top{}perTarget".format(n)]).union(context),
                transcription_factor=tf_name,
                gene2weights=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


ACTIVATING_MODULE = "activating"
REPRESSING_MODULE = "repressing"


def modules_from_adjacencies(adjacencies: pd.DataFrame,
                             ex_mtx: pd.DataFrame,
                        nomenclature: str,
                        thresholds=(0.001,0.005),
                        top_n_targets=(50,),
                        top_n_regulators=(5,10,50),
                        min_genes=20):
    """
    Create modules from a dataframe containing weighted adjacencies between a TF and a target genes.
    
    :param adjacencies: The dataframe with the TF-target links.
    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param nomenclature: The nomenclature of the genes.
    :param thresholds: the first method to create the TF-modules based on the best targets for each transcription factor.
    :param top_n_targets: the second method is to select the top targets for a given TF.
    :param top_n_regulators: the alternative way to create the TF-modules is to select the best regulators for each gene.
    :param min_genes: The required minimum number of genes in a module.
    :return: A list of Regulomes.
    """

    # Relationship between TF and its target, i.e. activator or repressor, is derived using the original expression
    # profiles. The Pearson product-moment correlation coefficient is used to derive this information.

    # Add correlation column and create two disjoint set of adjacencies.
    adjacencies = add_correlation(adjacencies.copy(), ex_mtx)
    activating_modules = adjacencies[adjacencies['correlation'] > 0.0]
    repressing_modules = adjacencies[adjacencies['correlation'] < 0.0]

    # Derive modules for these two sets of adjacencies
    # Add the transcription factor to the module.
    # Filter for minimum number of genes.
    def iter_modules(adjc, context):
        yield from chain(chain.from_iterable(modules4thr(adjc, thr, nomenclature, context) for thr in thresholds),
                         chain.from_iterable(modules4top_targets(adjc, n, nomenclature, context) for n in top_n_targets),
                         chain.from_iterable(modules4top_factors(adjc, n, nomenclature, context) for n in top_n_regulators))
    def add_tf(module):
        return module.add(module.transcription_factor)
    return list(filter(lambda m: len(m) >= min_genes,
                    map(add_tf,
                       chain(iter_modules(activating_modules, frozenset([ACTIVATING_MODULE])),
                             iter_modules(repressing_modules, frozenset([REPRESSING_MODULE]))))))


def save_to_yaml(signatures: Sequence[Type[GeneSignature]], fname: str):
    """

    :param signatures:
    :return:
    """
    with open(fname, 'w') as f:
        f.write(dump(signatures, default_flow_style=False, Dumper=Dumper))


def load_from_yaml(fname: str) -> Sequence[Type[GeneSignature]]:
    """

    :param fname:
    :return:
    """
    with open(fname, 'r') as f:
        return load(f.read(), Loader=Loader)


def add_motif_url(df: pd.DataFrame, base_url: str):
    """

    :param df:
    :param base_url:
    :return:
    """
    df[("Enrichment", COLUMN_NAME_MOTIF_URL)] = list(map(partial(urljoin, base=base_url), df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    return df

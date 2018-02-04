# -*- coding: utf-8 -*-

import pandas as pd
from .genesig import Regulome
from collections import defaultdict, Counter
from itertools import chain


COLUMN_NAME_MOTIF_SIMILARITY_QVALUE = 'motif_similarity_qvalue'
COLUMN_NAME_ORTHOLOGOUS_IDENTITY = 'orthologous_identity'


def load_motif_annotations(fname: str,
                           column_names=('#motif_id', 'gene_name',
                                         COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
                                         'description')) -> pd.DataFrame:
    """
    Load motif annotations from a motif2TF snapshot.

    :param fname: the snapshot taken from motif2TF.
    :param column_names: the names of the columns in the snapshot to load.
    :return: A dataframe.
    """
    # Create a MultiIndex for the index combining unique gene name and motif ID. This should facilitate
    # later merging.
    return pd.read_csv(fname, sep='\t', index_col=[1,0], usecols=column_names)


COLUMN_NAME_TF = "TF"
COLUMN_NAME_TARGET = "target"
COLUMN_NAME_WEIGHT = "importance"


def modules4thr(adjacencies, threshold, nomenclature="MGI"):
    """

    :param adjacencies:
    :param threshold:
    :param nomenclature:
    :return:
    """
    for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
        module = target_genes[target_genes[COLUMN_NAME_WEIGHT] >= threshold]
        if len(module) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=("target weight >= {}".format(threshold)),
                transcription_factor=tf_name,
                gene2weights=list(zip(module[COLUMN_NAME_TARGET].values, module[COLUMN_NAME_WEIGHT].values)))


def modules4top_targets(adjacencies, n, nomenclature="MGI"):
    """

    :param adjacencies:
    :param n:
    :param nomenclature:
    :return:
    """
    for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
        module = target_genes.sort_values(by=COLUMN_NAME_WEIGHT, ascending=False).head(n)
        if len(module) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=("target in top {}".format(n)),
                transcription_factor=tf_name,
                gene2weights=list(zip(module[COLUMN_NAME_TARGET].values, module[COLUMN_NAME_WEIGHT].values)))


def modules4top_factors(adjacencies, n, nomenclature="MGI"):
    """

    :param adjacencies:
    :param n:
    :param nomenclature:
    :return:
    """
    tf2target2weight = defaultdict(Counter)
    for target_name, factors in adjacencies.groupby(by=COLUMN_NAME_TARGET):
        regulators = factors.sort_values(by=COLUMN_NAME_WEIGHT, ascending=False).head(n)
        for factor, weight in zip(regulators[COLUMN_NAME_TF].values, regulators[COLUMN_NAME_WEIGHT].values):
            tf2target2weight[factor][target_name] = weight

    for tf_name, target2weight in tf2target2weight.items():
        yield Regulome(
            name="Regulome for {}".format(tf_name),
            nomenclature=nomenclature,
            context=("factor in top {}".format(n)),
            transcription_factor=tf_name,
            gene2weights=target2weight)


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
# -*- coding: utf-8 -*-

from pandas import DataFrame
from .genesig import Regulome
from collections import defaultdict, Counter
from itertools import chain


COLUMN_NAME_TF = "TF"
COLUMN_NAME_TARGET = "target"
COLUMN_NAME_WEIGHT = "importance"


def regulome4thr(adjacencies, threshold, nomenclature="MGI"):
    """

    :param adjacencies:
    :param threshold:
    :param nomenclature:
    :return:
    """
    for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
        regulome = target_genes[target_genes[COLUMN_NAME_WEIGHT] >= threshold]
        if len(regulome) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=("target weight >= {}".format(threshold)),
                transcription_factor=tf_name,
                gene2weights=list(zip(regulome[COLUMN_NAME_TARGET].values, regulome[COLUMN_NAME_WEIGHT].values)))


def regulome4top_targets(adjacencies, n, nomenclature="MGI"):
    """

    :param adjacencies:
    :param n:
    :param nomenclature:
    :return:
    """
    for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
        regulome = target_genes.sort_values(by=COLUMN_NAME_WEIGHT, ascending=False).head(n)
        if len(regulome) > 0:
            yield Regulome(
                name="Regulome for {}".format(tf_name),
                nomenclature=nomenclature,
                context=("target in top {}".format(n)),
                transcription_factor=tf_name,
                gene2weights=list(zip(regulome[COLUMN_NAME_TARGET].values, regulome[COLUMN_NAME_WEIGHT].values)))


def regulome4top_factors(adjacencies, n, nomenclature="MGI"):
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


def regulomes_from_genie3(adjacencies: DataFrame, nomenclature: str,
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
    yield from chain(chain.from_iterable(regulome4thr(adjacencies, thr, nomenclature) for thr in thresholds),
                     chain.from_iterable(regulome4top_targets(adjacencies, n, nomenclature) for n in top_n_targets),
                     chain.from_iterable(regulome4top_factors(adjacencies, n, nomenclature) for n in top_n_regulators))
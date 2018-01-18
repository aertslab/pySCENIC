# -*- coding: utf-8 -*-

from pandas import DataFrame
from .genesig import Regulome


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
    COLUMN_NAME_TF = "TF"
    COLUMN_NAME_TARGET = "targets"
    COLUMN_NAME_WEIGHT = "importance"


    #TODO: Adjust to weights!

    def regulome4thr(threshold):
        for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
            return Regulome(
                name="Regulome for {} (target weight >= {})".format(tf_name, threshold),
                nomenclature=nomenclature,
                transcription_factor=tf_name,
                genes=target_genes[target_genes[COLUMN_NAME_WEIGHT] >= threshold][COLUMN_NAME_TARGET].values)

    def regulome4top_targets(n):
        for tf_name, target_genes in adjacencies.groupby(by=COLUMN_NAME_TF):
            return Regulome(
                name="Regulome for {} (target in top {})".format(tf_name, n),
                nomenclature=nomenclature,
                transcription_factor=tf_name,
                genes=target_genes.sort_values(by=COLUMN_NAME_WEIGHT, ascending=False).head(n)[COLUMN_NAME_TARGET].values)

    def regulome4top_factors(n):
        for target_name, factors in adjacencies.groupby(by=COLUMN_NAME_TARGET):


            return Regulome(
                name="Regulome for {} (factor in top {})".format(tf_name, n),
                nomenclature=nomenclature,
                transcription_factor=tf_name,
                genes=target_genes.sort_values(by=COLUMN_NAME_WEIGHT, ascending=False).head(n)[COLUMN_NAME_TARGET].values)

    yield None
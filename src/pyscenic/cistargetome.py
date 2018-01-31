# -*- coding: utf-8 -*-

from .recovery import enrichment, leading_edge
from .genesig import Regulome

import pandas as pd
from functools import partial
from itertools import repeat


def generate_features(db, gs, rank_threshold=1500):
    """
    Find the enriched features
    """
    df = enrichment(db, gs, rank_threshold=rank_threshold).sort_values(by=('Enrichment', 'NES'), ascending=False)
    df[('Metadata', 'Signature')] = gs.name
    df[('Metadata', 'Database')] = db.name
    if isinstance(gs, Regulome):
        df[('Metadata', 'Factor')] = gs.transcription_factor
    return df


def generate_recovery_curves(df_features):
    """

    """
    avgrcc = df_features['Recovery'].mean(axis=0)
    stdrcc = df_features['Recovery'].std(axis=0)
    return avgrcc + 2.0 * stdrcc


def filter_features(df_features, nes_threshold=3.0):
    """

    :param df_features:
    :param nes_threshold:
    :return:
    """
    return df_features[df_features[('Enrichment', 'NES')] >= nes_threshold]


def load_motif2tf_snapshot(fname):
    COLUMN_NAMES = ['gene_name', 'motif_similarity_qvalue', 'orthologous_identity', 'description']
    motif2tf = pd.read_csv(fname, sep='\t', index_col=0)
    motif2tf = motif2tf[COLUMN_NAMES]
    motif2tf.columns = list(zip(repeat('Motif2TF'), COLUMN_NAMES))
    return motif2tf


def add_tf_annotations(df_features, motif2tf):
    df = pd.merge(df_features, motif2tf, how='left', left_index=True, right_index=True)
    df = df[df[('Metadata', 'Factor')] == df['Motif2TF', 'gene_name']]
    return df


def add_targetome(df_features, avg2stdrcc, nomenclature="MGI"):
    df_features[('Enrichment', 'LE')] = df_features[['Recovery', 'Ranking']].apply(
        partial(leading_edge,
                avg2stdrcc=avg2stdrcc,
                genes=df_features['Ranking'].columns.values,
                nomenclature=nomenclature),
        axis=1)
    del df_features['Ranking']
    return df_features




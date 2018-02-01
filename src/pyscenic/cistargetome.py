# -*- coding: utf-8 -*-

from .recovery import enrichment, leading_edge
from .genesig import Regulome

import pandas as pd
from functools import partial
from itertools import repeat
from typing import Type
from .genesig import GeneSignature
from .rnkdb import RankingDatabase
import math
from operator import itemgetter
from itertools import chain
from dask.multiprocessing import get
from dask import delayed


def load_motif2tf_snapshot(fname: str) -> pd.DataFrame:
    """
    Load a motif2TF snapshot.
    """
    COLUMN_NAMES = ['gene_name', 'motif_similarity_qvalue', 'orthologous_identity', 'description']

    motif2tf = pd.read_csv(fname, sep='\t', index_col=0)
    motif2tf = motif2tf[COLUMN_NAMES]
    motif2tf.columns = pd.MultiIndex.from_tuples(list(map(tuple, zip(repeat('Motif2TF'), COLUMN_NAMES))))
    return motif2tf


def generate_features(db: RankingDatabase, gs: Type[GeneSignature], rank_threshold: int =1500) -> pd.DataFrame:
    """
    Find the enriched features from the supplied database for the given gene signature.
    """
    assert db
    assert gs
    assert rank_threshold >= 1

    # Create dataframe with enriched regulatory features.
    df = enrichment(db, gs, rank_threshold=rank_threshold).sort_values(by=('Enrichment', 'NES'), ascending=False)

    # Add additional metadata regarding the signature to the dataframe.
    df[('Metadata', 'Signature')] = gs.name
    df[('Metadata', 'Database')] = db.name
    df[('Metadata', 'Factor')] = gs.transcription_factor if isinstance(gs, Regulome) else None

    return df


def generate_recovery_curves(df_features: pd.DataFrame) -> pd.DataFrame:
    """
    Generate recovery curves
    """
    assert len(df_features) > 0

    avgrcc = df_features['Recovery'].mean(axis=0)
    stdrcc = df_features['Recovery'].std(axis=0)

    return pd.DataFrame(data={"avg": avgrcc, "avg2std": avgrcc + 2.0 * stdrcc})


def filter_features(df_features: pd.DataFrame, nes_threshold=3.0) -> pd.DataFrame:
    """
    Only keep significant regulatory features.
    """
    return df_features[df_features[('Enrichment', 'NES')] >= nes_threshold]


def add_tf_annotations(df_features: pd.DataFrame, motif2tf: pd.DataFrame):
    """
    Add TF annotations for enriched regulatory features.
    """
    return pd.merge(df_features, motif2tf, how='left', left_index=True, right_index=True)


def filter_annotations(df_features: pd.DataFrame) -> pd.DataFrame:
    """
    Filter enriched features on correct TF annotations.
    """
    return df_features[df_features[('Metadata', 'Factor')] == df_features[('Motif2TF', 'gene_name')]]


def add_targetome(df_features: pd.DataFrame, avg2stdrcc: pd.Series) -> pd.DataFrame:
    """
    Add targetome for enriched regulatory features.
    """
    df_features[('Enrichment', 'LE')] = df_features[['Recovery', 'Ranking']].apply(
        partial(leading_edge,
                avg2stdrcc=avg2stdrcc,
                genes=df_features['Ranking'].columns.values),
        axis=1)
    del df_features['Ranking']
    return df_features


def add_regulome_score(df_features: pd.DataFrame) -> pd.DataFrame:
    """
    Score the regulome taking into account NES of the motif and the strength of the TF annotation.
    """
    def score4row(row):
        nes = row[('Enrichment', 'NES')]
        motif_similarity_qval = row[('Motif2TF', 'motif_similarity_qvalue')]
        orthologuous_identity = row[('Motif2TF', 'orthologous_identity')]
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    df_features[('Enrichment', 'Score')] = df_features.apply(score4row, axis=1)
    return df_features


def create_regulomes(df_features, nomenclature="MGI"):
    def combine_leading_edges(les):
        return frozenset(chain.from_iterable(map(itemgetter(0), le) for le in les))

    def regulomes():
        for metadata, group in df_features.groupby(by=[('Metadata', 'Signature'),
                                                       ('Metadata', 'Database'),
                                                       ('Metadata', 'Factor')]):
            tf_name = metadata[2]
            regulome_name = "Regulome for {}".format(metadata[0], metadata[1])
            signature = combine_leading_edges(group[('Enrichment', 'LE')].values)
            name =  metadata[0]
            context = (metadata[1], name[name.find("(")+1:-2])
            score = group[("Enrichment", "Score")].max()
            yield Regulome(name=regulome_name,
                           transcription_factor=tf_name,
                           nomenclature=nomenclature,
                           score=score,
                           context=context,
                           gene2weights=signature)

    return list(regulomes())


def cistargetome(dbs, modules, fname, rank_threshold, nes_threshold, num_workers):
    """

    :param dbs:
    :param modules:
    :param fname:
    :param rank_threshold:
    :param nes_threshold:
    :param num_workers:
    :return:
    """
    dsk_generate_features = delayed(generate_features)
    dsk_generate_recovery_curves = delayed(generate_recovery_curves)
    dsk_add_targetome = delayed(add_targetome)

    @delayed
    def dsk_annotate_features(df_features, motif2tf, nes_threshold):
        return add_regulome_score(filter_annotations(add_tf_annotations(filter_features(df_features, nes_threshold), motif2tf)))

    @delayed
    def dsk_derive_regulomes(dfs):
        return create_regulomes(pd.concat(dfs))\

    motif2tf = load_motif2tf_snapshot(fname)
    features = [dsk_generate_features(db, gs, rank_threshold=rank_threshold) for db in dbs for gs in modules]
    rccs = [ dsk_generate_recovery_curves(f) for f in features ]
    annot_features = [dsk_annotate_features(f, motif2tf, nes_threshold) for f in features]
    targetomes = [dsk_add_targetome(af, rcc) for af, rcc in zip(annot_features, rccs) ]
    regulomes = dsk_derive_regulomes(targetomes)

    return regulomes.compute(get=get, num_workers=num_workers)

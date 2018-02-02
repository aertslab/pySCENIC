# -*- coding: utf-8 -*-

from .recovery import recovery, leading_edge
import pandas as pd
import numpy as np
from .utils import load_motif2tf_snapshot
from itertools import repeat
from .rnkdb import RankingDatabase
from functools import reduce, partial
from operator import concat
from dask.multiprocessing import get
from dask import delayed
import multiprocessing
from typing import Type, List, Sequence
from .genesig import Regulome
from .genesig import GeneSignature
from collections import defaultdict
import math


COLUMN_NAME_NES = "NES"
COLUMN_NAME_AUC = "AUC"


def derive_regulomes4pair(db: Type[RankingDatabase], regulome: Regulome, motif_annotations_fname: str,
                          rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0) -> List[Regulome]:
    """
    Calculate all regulomes for a given ranking database and a co-expression module.

    :param db: The ranking database.
    :param regulome: The co-expression module.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) for
    :param num_workers: The number of workers to use for the calculation. None of all available CPUs need to be used.
    :return: A list of regulomes.
    """
    #TODO: Pass database as filename? Will this result in a performance improvement?
    #TODO: Pass motif2TF as preloaded dataframe or load each time again? How can data be shared between processes?

    # Load rank of genes from database.
    df = db.load(regulome)
    features, genes, rankings = df.index.values, df.columns.values, df.values
    weights = np.asarray([regulome[gene] for gene in genes])

    # Calculate recovery curves, AUC and NES values.
    rccs, aucs = recovery(df, db.total_genes, weights, rank_threshold, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(features[enriched_features_idx],
                                                                              repeat(regulome.transcription_factor)))),
                                     data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                           COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    if len(enriched_features) == 0:
        return []

    # Find motif annotations for enriched features.
    #TODO: Potential performance improvement: shared memory version anongst processes (i.e. dask tasks). If sharing not
    #TODO: possible then penalty is pickling of this dataframe between processes.
    motif_annotations = load_motif2tf_snapshot(motif_annotations_fname)
    #TODO: Potential performance improvements:
    #TODO: 1. Invert levels of index MultiIndex?
    #TODO: 2. Reduce size of the motif annotations by filtering out weak annotations or only keep the best annotation
    #TODO:    for each gene.
    #TODO: 3. Write a custom dict-based implementation: (motif_ID, gene_name) => annotations OR
    #TODO:    gene_name => motif_ID => annotations [CAVE: dictionary creation time might be penalty!]
    annotated_features = pd.merge(enriched_features, motif_annotations, how="inner", left_index=True, right_index=True)
    if len(annotated_features) == 0:
        return []

    # Calculated leading edge for the remaining enriched features that have annotations.
    avgrcc = rccs.mean(axis=0)
    avg2stdrcc =  avgrcc + 2.0 * rccs.std(axis=0)
    context2regulomes = defaultdict(set)

    def score(nes, motif_similarity_qval, orthologuous_identity):
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) and motif_similarity_qval != 0.0 else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    for idx, row, rcc, ranking in zip(annotated_features.iterrows(), rccs[enriched_features_idx, :], rankings[enriched_features_idx, :]):
        context = (regulome.context[0], db.name)
        context2regulomes[context].add(Regulome(name=regulome.name,
                       score=score(row[COLUMN_NAME_NES], row['motif_similarity_qvalue'], row['orthologous_identity']),
                       nomenclature=regulome.nomenclature,
                       context=(regulome.context[0], db.name),
                       transcription_factor=regulome.transcription_factor,
                       gene2weights=leading_edge(rcc, avg2stdrcc, ranking, genes, regulome)))

    # Group regulomes per context.
    return [reduce(Regulome.union, regulomes) for regulomes in context2regulomes.values()]


def derive_regulomes(rnkdbs: Sequence[Type[RankingDatabase]], modules: Sequence[Type[GeneSignature]],
                         rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                         num_workers=None) -> List[Regulome]:
    """
    Calculate all regulomes for a given sequence of ranking databases and a sequence of co-expression modules.

    :param rnkdbs: The sequence of databases.
    :param modules: The sequence of modules.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) for
    :param num_workers: The number of workers to use for the calculation. None of all available CPUs need to be used.
    :return: A list of regulomes.
    """
    regulomes = delayed(partial(reduce, func=concat))([delayed(derive_regulomes4pair)(db, gs,
                                                                                 rank_threshold, auc_threshold, nes_threshold)
                                                       for db in rnkdbs for gs in modules])
    return regulomes.compute(get=get, num_workers=num_workers if num_workers else multiprocessing.cpu_count())

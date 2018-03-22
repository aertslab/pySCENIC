# -*- coding: utf-8 -*-

from .recovery import recovery, aucs as calc_aucs
import logging
import traceback
import pandas as pd
import numpy as np
from .utils import COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY, \
    COLUMN_NAME_MOTIF_ID, COLUMN_NAME_TF, COLUMN_NAME_ANNOTATION, ACTIVATING_MODULE, REPRESSING_MODULE
from itertools import repeat
from .rnkdb import RankingDatabase
from functools import reduce
from typing import Type, Sequence, Optional
from .genesig import Regulome
from .recovery import leading_edge4row
import math
from itertools import chain
from functools import partial
from cytoolz import first
import numpy as np
from dask.dataframe.utils import make_meta


COLUMN_NAME_NES = "NES"
COLUMN_NAME_AUC = "AUC"
COLUMN_NAME_CONTEXT = "Context"
COLUMN_NAME_TARGET_GENES = "TargetGenes"
COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
COLUMN_NAME_TYPE = "Type"
#TODO: Should actually be a function depending on return_recovery_curves and rank_threshold
DF_META_DATA = make_meta({('Enrichment', COLUMN_NAME_AUC): np.float64,
                          ('Enrichment', COLUMN_NAME_NES): np.float64,
                          ('Enrichment', COLUMN_NAME_MOTIF_SIMILARITY_QVALUE): np.float64,
                          ('Enrichment', COLUMN_NAME_ORTHOLOGOUS_IDENTITY): np.float64,
                          ('Enrichment', COLUMN_NAME_ANNOTATION): np.object,
                          ('Enrichment', COLUMN_NAME_CONTEXT): np.object,
                          ('Enrichment', COLUMN_NAME_TARGET_GENES): np.object,
                          ('Enrichment', COLUMN_NAME_RANK_AT_MAX): np.int64})


__all__ = ["module2features", "module2df", "modules2df", "df2regulomes", "module2regulome", "modules2regulomes"]


LOGGER = logging.getLogger(__name__)


def module2features_rcc4all_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                                 rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                                 avgrcc_sample_frac: float = None, weighted_recovery=False,
                                 filter_for_annotation=True):
    """
    Create a dataframe of enriched and annotated features a given ranking database and a co-expression module.

    :param db: The ranking database.
    :param module: The co-expression module.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) threshold to select enriched features.
    :param avgrcc_sample_frac: The fraction of the features to use for the calculation of the average curve, If None
        then all features are used.
    :param weighted_recovery: Use weighted recovery in the analysis.
    :return: A dataframe with enriched and annotated features.
    """

    # Load rank of genes from database.
    df = db.load(module)
    features, genes, rankings = df.index.values, df.columns.values, df.values
    weights = np.asarray([module[gene] for gene in genes]) if weighted_recovery else np.ones(len(genes))

    # Calculate recovery curves, AUC and NES values.
    rccs, aucs = recovery(df, db.total_genes, weights, rank_threshold, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(repeat(module.transcription_factor),
                                                                              features[enriched_features_idx])),
                                                                     names=[COLUMN_NAME_TF, COLUMN_NAME_MOTIF_ID]),
                                     data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                           COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    if len(enriched_features) == 0:
        return pd.DataFrame(), None, None, genes, None

    # Find motif annotations for enriched features.
    annotated_features = pd.merge(enriched_features, motif_annotations, how="left", left_index=True, right_index=True)
    annotated_features_idx = pd.notnull(annotated_features[COLUMN_NAME_ANNOTATION]) if filter_for_annotation else np.full((len(enriched_features),), True)
    if len(annotated_features[annotated_features_idx]) == 0:
        return pd.DataFrame(), None, None, genes, None

    # Calculated leading edge for the remaining enriched features that have annotations.
    if avgrcc_sample_frac is None:
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc =  avgrcc + 2.0 * rccs.std(axis=0)
    else:
        n_features = len(features)
        sample_idx = np.random.randint(0, int(n_features*avgrcc_sample_frac))
        avgrcc = rccs[sample_idx, :].mean(axis=0)
        avg2stdrcc = avgrcc + 2.0 * rccs[sample_idx, :].std(axis=0)
    rccs = rccs[enriched_features_idx, :][annotated_features_idx, :]
    rankings = rankings[enriched_features_idx, :][annotated_features_idx, :]

    # Add additional information to the dataframe.
    annotated_features = annotated_features[annotated_features_idx]
    context = frozenset(chain(module.context, [db.name]))
    annotated_features[COLUMN_NAME_CONTEXT] = len(annotated_features) * [context]

    return annotated_features, rccs, rankings, genes, avg2stdrcc


def module2features_auc1st_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                                rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                                avgrcc_sample_frac: float = None, weighted_recovery=False,
                                filter_for_annotation=True):
    """
    Create a dataframe of enriched and annotated features a given ranking database and a co-expression module.

    :param db: The ranking database.
    :param module: The co-expression module.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) threshold to select enriched features.
    :param avgrcc_sample_frac: The fraction of the features to use for the calculation of the average curve, If None
        then all features are used.
    :param weighted_recovery: Use weighted recovery in the analysis.
    :return: A dataframe with enriched and annotated features.
    """

    # Load rank of genes from database.
    df = db.load(module)
    features, genes, rankings = df.index.values, df.columns.values, df.values
    weights = np.asarray([module[gene] for gene in genes]) if weighted_recovery else np.ones(len(genes))

    # Calculate recovery curves, AUC and NES values.
    # For fast unweighted implementation so weights to None.
    aucs = calc_aucs(df, db.total_genes, weights, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(repeat(module.transcription_factor),
                                                                              features[enriched_features_idx])),
                                                                     names=[COLUMN_NAME_TF, COLUMN_NAME_MOTIF_ID]),
                                     data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                           COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    if len(enriched_features) == 0:
        return pd.DataFrame(), None, None, genes, None

    # Find motif annotations for enriched features.
    annotated_features = pd.merge(enriched_features, motif_annotations, how="left", left_index=True, right_index=True)
    annotated_features_idx = pd.notnull(annotated_features[COLUMN_NAME_ANNOTATION]) if filter_for_annotation else np.full((len(enriched_features),), True)
    if len(annotated_features[annotated_features_idx]) == 0:
        return pd.DataFrame(), None, None, genes, None

    # Calculated leading edge for the remaining enriched features that have annotations.
    if avgrcc_sample_frac is None:
        rccs, _ = recovery(df, db.total_genes, np.full(len(genes), 1.0), rank_threshold, auc_threshold, no_auc=True)
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc =  avgrcc + 2.0 * rccs.std(axis=0)
    else:
        rccs, _ = recovery(df.sample(frac=avgrcc_sample_frac), db.total_genes, np.full(len(genes), 1.0), rank_threshold, auc_threshold, no_auc=True)
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc = avgrcc + 2.0 * rccs.std(axis=0)
    rccs = rccs[enriched_features_idx, :][annotated_features_idx, :]
    rankings = rankings[enriched_features_idx, :][annotated_features_idx, :]

    # Add additional information to the dataframe.
    annotated_features = annotated_features[annotated_features_idx]
    context = frozenset(chain(module.context, [db.name]))
    annotated_features[COLUMN_NAME_CONTEXT] = len(annotated_features) * [context]

    return annotated_features, rccs, rankings, genes, avg2stdrcc


module2features = partial(module2features_auc1st_impl,
                          rank_threshold = 1500, auc_threshold = 0.05, nes_threshold=3.0,
                          avgrcc_sample_frac = None, filter_for_annotation=True)


def module2df(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
              weighted_recovery=False, return_recovery_curves=False, module2features_func=module2features) -> pd.DataFrame:
    """

    """
    # Derive enriched and TF-annotated features for module.
    try:
        df_annotated_features, rccs, rankings, genes, avg2stdrcc = module2features_func(db, module, motif_annotations,
                                                                                    weighted_recovery=weighted_recovery)
    except MemoryError:
        LOGGER.error("Unable to process \"{}\" on database \"{}\" because ran out of memory. Stacktrace:".format(module.name, db.name))
        LOGGER.error(traceback.format_exc())
        return DF_META_DATA
    # If less than 80% of the genes are mapped to the ranking database, the module is skipped.
    n_missing = len(module) - len(genes)
    frac_missing = float(n_missing)/len(module)
    if frac_missing >= 0.20:
        LOGGER.warning("Less than 80% of the genes in {} could be mapped to {}. Skipping this module.".format(module.name, db.name))
        return DF_META_DATA

    # If no annotated enriched features could be found, skip module.
    if len(df_annotated_features) == 0:
        return DF_META_DATA
    rank_threshold = rccs.shape[1]

    # Combine elements into a dataframe.
    df_annotated_features.columns = pd.MultiIndex.from_tuples(list(zip(repeat("Enrichment"),
                                                                       df_annotated_features.columns)))
    df_rnks = pd.DataFrame(index=df_annotated_features.index,
                           columns=pd.MultiIndex.from_tuples(list(zip(repeat("Ranking"), genes))),
                           data=rankings)
    df_rccs = pd.DataFrame(index=df_annotated_features.index,
                           columns=pd.MultiIndex.from_tuples(list(zip(repeat("Recovery"), np.arange(rank_threshold)))),
                           data=rccs)
    df = pd.concat([df_annotated_features, df_rccs, df_rnks], axis=1)

    # Calculate the leading edges for each row. Always return importance from gene inference phase.
    weights = np.array([module[gene] for gene in genes])
    df[[("Enrichment", COLUMN_NAME_TARGET_GENES), ("Enrichment", COLUMN_NAME_RANK_AT_MAX)]] = df.apply(partial(leading_edge4row,
                                                                                                               avg2stdrcc=avg2stdrcc, genes=genes, weights=weights), axis=1)

    # Remove unnecessary data from dataframe.
    del df['Ranking']
    if not return_recovery_curves:
        del df['Recovery']
    return df


def modules2df(db: Type[RankingDatabase], modules: Sequence[Regulome], motif_annotations: pd.DataFrame,
               weighted_recovery=False, module2features_func=module2features) -> pd.DataFrame:
    # Make sure return recovery curves is always set to false because the metadata for the distributed dataframe needs
    # to be fixed for the dask framework.
    #TODO: Remove this restriction.
    return pd.concat([module2df(db, module, motif_annotations, weighted_recovery, False, module2features_func)
                      for module in modules])


def regulome4group(tf_name, context, df_group, nomenclature) -> Optional[Regulome]:
    """

    """
    # Create regulomes for each enriched and annotated feature.
    def score(nes, motif_similarity_qval, orthologuous_identity):
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) and motif_similarity_qval != 0.0 else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    def derive_interaction_type(ctx):
        return "(-)" if REPRESSING_MODULE in ctx else "(+)"

    def row2regulome(row):
        # The target genes as well as their weights/importances are directly taken from the dataframe.
        return Regulome(name="{}{}".format(tf_name,derive_interaction_type(context)),
                        score=score(row[COLUMN_NAME_NES],
                                    row[COLUMN_NAME_MOTIF_SIMILARITY_QVALUE],
                                    row[COLUMN_NAME_ORTHOLOGOUS_IDENTITY]),
                        nomenclature=nomenclature,
                        context=context,
                        transcription_factor=tf_name,
                        gene2weights=row[COLUMN_NAME_TARGET_GENES])

    # Aggregate these regulomes into a single one using the union operator.
    return reduce(Regulome.union, (row2regulome(row) for _, row in df_group.iterrows()))


def df2regulomes(df, nomenclature) -> Sequence[Regulome]:
    """
    Create regulomes from a dataframe of enriched features.

    :param df:
    :param nomenclature: The nomenclature of the genes.
    :return: A sequence of regulomes.
    """
    # Normally the columns index has two levels. For convenience of the following code the first level is removed.
    if df.columns.nlevels == 2:
        df.columns = df.columns.droplevel(0)

    def get_type(row):
        ctx = row[COLUMN_NAME_CONTEXT]
        return ACTIVATING_MODULE if ACTIVATING_MODULE in ctx else REPRESSING_MODULE
    df[COLUMN_NAME_TYPE] = df.apply(get_type,axis=1)

    not_none = lambda r: r is not None
    return list(filter(not_none, (regulome4group(tf_name, frozenset([interaction_type]), df_grp, nomenclature)
                                        for (tf_name, interaction_type), df_grp in df.groupby(by=[COLUMN_NAME_TF,
                                                                                                  COLUMN_NAME_TYPE]))))


def module2regulome(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                    weighted_recovery=False, return_recovery_curves=False,
                    module2features_func=module2features) -> Optional[Regulome]:
    # First calculating a dataframe and then derive the regulomes from them introduces a performance penalty.
    df = module2df(db, module, motif_annotations, weighted_recovery=weighted_recovery,
                   return_recovery_curves=return_recovery_curves,
                   module2features_func=module2features_func)
    if len(df) == 0:
        return None
    regulomes = df2regulomes(df, module.nomenclature)
    return first(regulomes) if len(regulomes) > 0 else None


def modules2regulomes(db: Type[RankingDatabase], modules: Sequence[Regulome], motif_annotations: pd.DataFrame,
                      weighted_recovery=False, return_recovery_curves=False,
                      module2features_func=module2features) -> Sequence[Regulome]:
    assert len(modules) > 0

    nomenclature = modules[0].nomenclature
    df = modules2df(db, modules, motif_annotations, weighted_recovery=weighted_recovery,
                    return_recovery_curves=return_recovery_curves,
                    module2features_func=module2features_func)
    return [] if len(df) == 0 else df2regulomes(df, nomenclature)

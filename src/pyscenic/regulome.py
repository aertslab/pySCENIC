# -*- coding: utf-8 -*-

from .recovery import recovery, aucs as calc_aucs
import pandas as pd
import numpy as np
from .utils import load_motif_annotations, COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY, COLUMN_NAME_MOTIF_ID, COLUMN_NAME_TF, COLUMN_NAME_ANNOTATION
from itertools import repeat
from .rnkdb import RankingDatabase, MemoryDecorator
from functools import reduce
from operator import concat
from boltons.iterutils import chunked_iter
from dask.multiprocessing import get
from dask import delayed
from dask.distributed import LocalCluster, Client
from typing import Type, Sequence, Optional, Union
from .genesig import Regulome, GeneSignature
from .recovery import leading_edge4row
import math
from itertools import chain
from math import ceil
from functools import partial
from cytoolz import first
# Using multiprocessing using dill package for pickling to avoid strange bugs.
from multiprocessing import cpu_count
from multiprocessing_on_dill.connection import Pipe
from multiprocessing_on_dill.context import Process
import datetime
from .utils import add_motif_url


COLUMN_NAME_NES = "NES"
COLUMN_NAME_AUC = "AUC"
COLUMN_NAME_CONTEXT = "Context"
COLUMN_NAME_TARGET_GENES = "TargetGenes"
COLUMN_NAME_RANK_AT_MAX = "RankAtMax"


__all__ = ["module2features", "module2df", "modules2df", "df2regulomes", "module2regulome", "modules2regulomes",
           "prune_targets", "find_motifs"]


def module2features_bincount_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
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


def module2features_numba_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
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
    aucs = calc_aucs(df, db.total_genes, weights if weighted_recovery else None, rank_threshold, auc_threshold)
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


module2features = partial(module2features_numba_impl,
                          rank_threshold = 1500, auc_threshold = 0.05, nes_threshold=3.0,
                          avgrcc_sample_frac = None, filter_for_annotation=True)


def module2df(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
              weighted_recovery=False, return_recovery_curves=False, module2features_func=module2features) -> pd.DataFrame:
    """

    """
    # Derive enriched and TF-annotated features for module.
    df_annotated_features, rccs, rankings, genes, avg2stdrcc = module2features_func(db, module, motif_annotations,
                                                                                    weighted_recovery=weighted_recovery)
    # If less than 80% of the genes are mapped to the ranking database, the module is skipped.
    n_missing = len(module) - len(genes)
    frac_missing = float(n_missing)/len(module)
    if frac_missing >= 0.20:
        print("Less than 80% of the genes in {} could be mapped to {}. Skipping this module.".format(module.name, db.name))
        return pd.DataFrame()

    # If no annotated enriched features could be found, skip module.
    if len(df_annotated_features) == 0:
        return pd.DataFrame()
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

    # Calculate the leading edges for each row. Rank is discarded.
    weights = np.array([module[gene] for gene in genes]) if weighted_recovery else np.ones(len(genes))
    df[[("Enrichment", COLUMN_NAME_TARGET_GENES), ("Enrichment", COLUMN_NAME_RANK_AT_MAX)]] = df.apply(partial(leading_edge4row,
                                                          avg2stdrcc=avg2stdrcc, genes=genes, weights=weights), axis=1)

    # Remove unnecessary data from dataframe.
    del df['Ranking']
    if not return_recovery_curves:
        del df['Recovery']
    return df


def modules2df(db: Type[RankingDatabase], modules: Sequence[Regulome], motif_annotations: pd.DataFrame,
               weighted_recovery=False, return_recovery_curves=False, module2features_func=module2features) -> pd.DataFrame:
    """

    """
    return pd.concat([module2df(db, module, motif_annotations, weighted_recovery, return_recovery_curves, module2features_func)
                      for module in modules])


def regulome4group(tf_name, context, df_group, nomenclature) -> Optional[Regulome]:
    """

    """
    # Create regulomes for each enriched and annotated feature.
    def score(nes, motif_similarity_qval, orthologuous_identity):
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) and motif_similarity_qval != 0.0 else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    def row2regulome(row):
        # The target genes as well as their weights/importances are directly taken from the dataframe.
        return Regulome(name=tf_name,
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

    """
    not_none = lambda r: r is not None
    return list(filter(not_none, (regulome4group(tf_name, frozenset(), df_grp['Enrichment'], nomenclature)
            for tf_name, df_grp in df.groupby(by=COLUMN_NAME_TF))))


def module2regulome(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                    weighted_recovery=False, return_recovery_curves=False,
                    module2features_func=module2features) -> Optional[Regulome]:
    """

    """
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
    """

    """
    assert len(modules) > 0

    nomenclature = modules[0].nomenclature
    df = modules2df(db, modules, motif_annotations, weighted_recovery=weighted_recovery,
                   return_recovery_curves=return_recovery_curves,
                   module2features_func=module2features_func)
    return [] if len(df) == 0 else df2regulomes(df, nomenclature)


def _prepare_client(client_or_address):
    """
    :param client_or_address: one of:
           * None
           * verbatim: 'local'
           * string address
           * a Client instance
    :return: a tuple: (Client instance, shutdown callback function).
    :raises: ValueError if no valid client input was provided.
    """
    # Credits to Thomas Moerman (arboretum package):
    # https://github.com/tmoerman/arboretum/blob/b065c6eade325ace104b2bb772ad15c78d573b1b/arboretum/algo.py#L139-L185

    if client_or_address is None or str(client_or_address).lower() == 'local':
        local_cluster = LocalCluster(diagnostics_port=None)
        client = Client(local_cluster)

        def close_client_and_local_cluster(verbose=False):
            if verbose:
                print('shutting down client and local cluster')

            client.close()
            local_cluster.close()

        return client, close_client_and_local_cluster

    elif isinstance(client_or_address, str) and client_or_address.lower() != 'local':
        client = Client(client_or_address)

        def close_client(verbose=False):
            if verbose:
                print('shutting down client')

            client.close()

        return client, close_client

    elif isinstance(client_or_address, Client):

        def close_dummy(verbose=False):
            if verbose:
                print('not shutting down client, client was created externally')

            return None

        return client_or_address, close_dummy

    else:
        raise ValueError("Invalid client specified {}".format(str(client_or_address)))


class Worker(Process):
    def __init__(self, name: str, db: Type[RankingDatabase], modules: Sequence[Regulome],
                 motif_annotations_fname: str, sender,
                 motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0,
                 transformation_func=module2regulome):
        super().__init__(name=name)
        self.database = db
        self.modules = modules
        self.motif_annotations_fname = motif_annotations_fname
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologuous_identity_threshold = orthologuous_identity_threshold
        self.transform_fnc = transformation_func
        self.sender = sender

    def run(self):
        # Load ranking database in memory.
        rnkdb = MemoryDecorator(self.database)
        print("{} - Worker {}: database loaded in memory.".format(datetime.datetime.now(), self.name))

        # Load motif annotations in memory.
        motif_annotations = load_motif_annotations(self.motif_annotations_fname,
                                                   motif_similarity_fdr=self.motif_similarity_fdr,
                                                   orthologous_identity_threshold=self.orthologuous_identity_threshold)
        print("{} - Worker {}: motif annotations loaded in memory.".format(datetime.datetime.now(), self.name))

        # Apply transformation on all modules.
        output = self.transform_fnc(rnkdb, self.modules, motif_annotations=motif_annotations)
        print("{} - Worker {}: All regulomes derived.".format(datetime.datetime.now(), self.name))

        # Sending information back to parent process.
        # Another approach might be to write a CSV file (for dataframes) or YAML file (for regulomes) to a temp.
        # file and share the name of the file with the parent process.
        # Serialization introduces a performance penalty!
        self.sender.send(output)
        self.sender.close()
        print("{} - Worker {}: Done.".format(datetime.datetime.now(), self.name))


def _distributed_calc(rnkdbs: Sequence[Type[RankingDatabase]], modules: Sequence[Type[GeneSignature]],
                      motif_annotations_fname: str,
                      transform_func, aggregate_func,
                      motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0,
                      client_or_address='custom_multiprocessing',
                      num_workers=None, module_chunksize=100) -> Union[Sequence[Regulome], pd.DataFrame]:
    def is_valid(client_or_address):
        if isinstance(client_or_address, int):
            return True
        elif isinstance(client_or_address, str) and client_or_address in {"custom_multiprocessing", "dask_multiprocessing", "local"}:
            return True
        elif isinstance(client_or_address, Client):
            return True
        return False
    assert is_valid(client_or_address), "\"{}\"is not valid for parameter client_or_address.".format(client_or_address)

    if client_or_address == 'custom_multiprocessing':
        # This implementation overcomes the I/O-bounded performance by the dask-based parallelized/distributed version.
        # Each worker (subprocess) loads a dedicated ranking database and motif annotation table into its own memory
        # space before consuming module. The implementation of each worker uses the AUC-first numba JIT based implementation
        # of the algorithm.
        assert len(rnkdbs) <= num_workers if num_workers else cpu_count(), "The number of databases is larger than the number of cores."
        amplifier = int((num_workers if num_workers else cpu_count())/len(rnkdbs))
        print("Using {} workers.".format(len(rnkdbs) * amplifier))
        receivers = []
        for db in rnkdbs:
            for idx, chunk in enumerate(chunked_iter(modules, ceil(len(modules)/float(amplifier)))):
                sender, receiver = Pipe()
                receivers.append(receiver)
                Worker("{}({})".format(db.name, idx+1), db, chunk, motif_annotations_fname, sender,
                       motif_similarity_fdr, orthologuous_identity_threshold, transform_func).start()
        return aggregate_func([recv.recv() for recv in receivers])
    else:
        # Load motif annotations.
        motif_annotations = load_motif_annotations(motif_annotations_fname,
                                                   motif_similarity_fdr=motif_similarity_fdr,
                                                   orthologous_identity_threshold=orthologuous_identity_threshold)

        # Create dask graph.
        # For performance reasons we analyze multiple modules for a database in a single node of the dask graph:
        # One combination of a rankings database with a single gene signature takes approximate 1 sec to complete
        # in an isolated environment (cf. notebooks). Using dask with a dask graph where one node/task in the graph
        # corresponds to processing one single database and a single signature results in dire performance. When
        # resorting to a dask graph where a single node/task corresponds to processing multiple signatures on a database
        # greatly boost performance: emperically we noticed a drop in overall duration from +4h:30m to +30m.
        # Current (unvalidated) explanation: the overhead of the scheduler for a 1 second task is too high. But
        # this is refuted by this blog: http://matthewrocklin.com/blog/work/2016/05/05/performant-task-scheduling (overhead
        # of sceduler is neglectable even if tasks/nodes take only 1s to complete.
        dask_graph = delayed(aggregate_func)(
            (delayed(transform_func)
             (db, gs_chunk, motif_annotations) for db in rnkdbs for gs_chunk in chunked_iter(modules, module_chunksize)))

        # Compute dask graph.
        if client_or_address == "dask_multiprocessing":
            return dask_graph.compute(get=get, num_workers=num_workers if num_workers else cpu_count())
        else:
            # Run via dask.distributed framework.
            client, shutdown_callback = _prepare_client(client_or_address)
            try:
                return client.compute(dask_graph)
            finally:
                shutdown_callback(False)


def find_motifs(rnkdbs: Sequence[Type[RankingDatabase]], signatures: Sequence[Type[GeneSignature]],
                motif_annotations_fname: str,
                rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0,
                avgrcc_sample_frac: float = None,
                weighted_recovery=False, client_or_address='custom_multiprocessing',
                num_workers=None, module_chunksize=100,
                motif_base_url: str = "http://motifcollections.aertslab.org/v9/") -> pd.DataFrame:
    module2features_func = partial(module2features_numba_impl,
                                   rank_threshold=rank_threshold,
                                   auc_threshold=auc_threshold,
                                   nes_threshold=nes_threshold,
                                   avgrcc_sample_frac=avgrcc_sample_frac,
                                   filter_for_annotation=False)
    transformation_func = partial(modules2df, module2features_func=module2features_func, weighted_recovery=weighted_recovery)
    aggregation_func = pd.concat
    df = _distributed_calc(rnkdbs, signatures, motif_annotations_fname, transformation_func, aggregation_func,
                      motif_similarity_fdr, orthologuous_identity_threshold, client_or_address,
                      num_workers, module_chunksize)
    return add_motif_url(df, base_url=motif_base_url)


def prune_targets(rnkdbs: Sequence[Type[RankingDatabase]], modules: Sequence[Regulome],
                         motif_annotations_fname: str,
                         rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                         motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0,
                         avgrcc_sample_frac: float = None,
                         weighted_recovery=False, client_or_address='custom_multiprocessing',
                         num_workers=None, module_chunksize=100, output="regulomes") -> Union[Sequence[Regulome], pd.DataFrame]:
    """
    Calculate all regulomes for a given sequence of ranking databases and a sequence of co-expression modules.

    :param rnkdbs: The sequence of databases.
    :param modules: The sequence of modules.
    :param motif_annotations_fname: The name of the file that contains the motif annotations to use.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) threshold to select enriched features.
    :param motif_similarity_fdr: The maximum False Discovery Rate to find factor annotations for enriched motifs.
    :param orthologuous_identity_threshold: The minimum orthologuous identity to find factor annotations
        for enriched motifs.
    :param avgrcc_sample_frac: The fraction of the features to use for the calculation of the average curve, If None
        then all features are used.
    :param weighted_recovery: Use weights of a gene signature when calculating recovery curves?
    :param num_workers: If not using a cluster, the number of workers to use for the calculation.
        None of all available CPUs need to be used.
    :param module_chunksize: The size of the chunk to use when using the dask framework.
    :param client_or_address: The client of IP address of the scheduler when working with dask. For local multi-core
        systems 'custom_multiprocessing' or 'dask_multiprocessing' can be supplied.
    :param output: The type of output requested, i.e. regulomes or a dataframe.
    :return: A sequence of regulomes.
    """
    assert output in {"regulomes", "df"}, "Invalid output type."

    module2features_func = partial(module2features_numba_impl,
                              rank_threshold=rank_threshold,
                              auc_threshold=auc_threshold,
                              nes_threshold=nes_threshold,
                              avgrcc_sample_frac=avgrcc_sample_frac,
                              filter_for_annotation=True)
    transformation_func = partial(modules2regulomes, module2features_func=module2features_func, weighted_recovery=weighted_recovery) \
        if output == "regulomes" else partial(modules2df, module2features_func=module2features_func, weighted_recovery=weighted_recovery)
    from toolz.curried import reduce
    aggregation_func = reduce(concat) if output == "regulomes" else pd.concat
    return _distributed_calc(rnkdbs, modules, motif_annotations_fname, transformation_func, aggregation_func,
                      motif_similarity_fdr, orthologuous_identity_threshold, client_or_address,
                      num_workers, module_chunksize)

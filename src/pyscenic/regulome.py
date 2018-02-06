# -*- coding: utf-8 -*-

from .recovery import recovery, leading_edge, aucs as calc_aucs
import pandas as pd
import numpy as np
from .utils import load_motif_annotations, COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY
from itertools import repeat
from .rnkdb import RankingDatabase, MemoryDecorator
from functools import reduce
from operator import concat
from dask.multiprocessing import get
from dask import delayed
from dask.distributed import LocalCluster, Client
from multiprocessing import Pipe, cpu_count, Process
from typing import Type, Sequence, Optional
from .genesig import Regulome
import math
from cytoolz import compose


COLUMN_NAME_NES = "NES"
COLUMN_NAME_AUC = "AUC"


def module2regulome_bincount_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                    rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                    avgrcc_sample_frac: float = None, weighted_recovery=False) -> Optional[Regulome]:
    """
    Create a regulome for a given ranking database and a co-expression module. If non can be created NoN is returned.

    :param db: The ranking database.
    :param module: The co-expression module.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) threshold to select enriched features.
    :param avgrcc_sample_frac: The fraction of the features to use for the calculation of the average curve, If None
        then all features are used.
    :param weighted_recovery: Use weights of a gene signature when calculating recovery curves?
    :return: A regulome or None.
    """

    # Load rank of genes from database.
    df = db.load(module)
    features, genes, rankings = df.index.values, df.columns.values, df.values
    weights = np.asarray([module[gene] for gene in genes])

    # Calculate recovery curves, AUC and NES values.
    rccs, aucs = recovery(df, db.total_genes, weights, rank_threshold, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(repeat(module.transcription_factor),
                                                                              features[enriched_features_idx])),
                                                                     names=["gene_name", "#motif_id"]),
                                     data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                           COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    if len(enriched_features) == 0:
        return None

    # Find motif annotations for enriched features.
    annotated_features = pd.merge(enriched_features, motif_annotations, how="inner", left_index=True, right_index=True)
    if len(annotated_features) == 0:
        return None

    # Calculated leading edge for the remaining enriched features that have annotations.
    if avgrcc_sample_frac is None:
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc =  avgrcc + 2.0 * rccs.std(axis=0)
    else:
        n_features = len(features)
        sample_idx = np.random.randint(0, int(n_features*avgrcc_sample_frac))
        avgrcc = rccs[sample_idx, :].mean(axis=0)
        avg2stdrcc = avgrcc + 2.0 * rccs[sample_idx, :].std(axis=0)

    # Create regulomes for each enriched and annotated feature.
    def score(nes, motif_similarity_qval, orthologuous_identity):
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) and motif_similarity_qval != 0.0 else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    regulomes = []
    _module = module if weighted_recovery else module.noweights()
    for (_, row), rcc, ranking in zip(annotated_features.iterrows(), rccs[enriched_features_idx, :], rankings[enriched_features_idx, :]):
        regulomes.append(Regulome(name=module.name,
                                  score=score(row[COLUMN_NAME_NES],
                                              row[COLUMN_NAME_MOTIF_SIMILARITY_QVALUE],
                                              row[COLUMN_NAME_ORTHOLOGOUS_IDENTITY]),
                                  nomenclature=module.nomenclature,
                                  context=module.context.union(frozenset([db.name])),
                                  transcription_factor=module.transcription_factor,
                                  gene2weights=leading_edge(rcc, avg2stdrcc, ranking, genes, _module)))

    # Aggregate these regulomes into a single one using the union operator.
    return reduce(Regulome.union, regulomes)


def module2regulome_numba_impl(db: Type[RankingDatabase], module: Regulome, motif_annotations: pd.DataFrame,
                    rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                    avgrcc_sample_frac: float = None) -> Optional[Regulome]:
    """
    Create a regulome for a given ranking database and a co-expression module. If non can be created NoN is returned.

    :param db: The ranking database.
    :param module: The co-expression module.
    :param rank_threshold: The total number of ranked genes to take into account when creating a recovery curve.
    :param auc_threshold: The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    :param nes_threshold: The Normalized Enrichment Score (NES) threshold to select enriched features.
    :param avgrcc_sample_frac: The fraction of the features to use for the calculation of the average curve, If None
        then all features are used.
    :return: A regulome or None.
    """

    # Load rank of genes from database.
    df = db.load(module)
    features, genes, rankings = df.index.values, df.columns.values, df.values

    # Calculate recovery curves, AUC and NES values.
    aucs = calc_aucs(df, db.total_genes, rank_threshold, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold
    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(repeat(module.transcription_factor),
                                                                              features[enriched_features_idx])),
                                                                     names=["gene_name", "#motif_id"]),
                                     data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                           COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    if len(enriched_features) == 0:
        return None

    # Find motif annotations for enriched features.
    annotated_features = pd.merge(enriched_features, motif_annotations, how="inner", left_index=True, right_index=True)
    if len(annotated_features) == 0:
        return None

    # Calculated leading edge for the remaining enriched features that have annotations.
    if avgrcc_sample_frac is None:
        rccs, _ = recovery(df, db.total_genes, np.full(len(genes), 1.0), rank_threshold, auc_threshold)
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc =  avgrcc + 2.0 * rccs.std(axis=0)
    else:
        rccs, _ = recovery(df.sample(frac=avgrcc_sample_frac), db.total_genes, np.full(len(genes), 1.0), rank_threshold, auc_threshold)
        avgrcc = rccs.mean(axis=0)
        avg2stdrcc = avgrcc + 2.0 * rccs.std(axis=0)

    # Create regulomes for each enriched and annotated feature.
    def score(nes, motif_similarity_qval, orthologuous_identity):
        MAX_VALUE = 100
        score = nes * -math.log(motif_similarity_qval)/MAX_VALUE if not math.isnan(motif_similarity_qval) and motif_similarity_qval != 0.0 else nes
        return score if math.isnan(orthologuous_identity) else score * orthologuous_identity

    regulomes = []
    for (_, row), rcc, ranking in zip(annotated_features.iterrows(), rccs[enriched_features_idx, :], rankings[enriched_features_idx, :]):
        regulomes.append(Regulome(name=module.name,
                                  score=score(row[COLUMN_NAME_NES],
                                              row[COLUMN_NAME_MOTIF_SIMILARITY_QVALUE],
                                              row[COLUMN_NAME_ORTHOLOGOUS_IDENTITY]),
                                  nomenclature=module.nomenclature,
                                  context=module.context.union(frozenset([db.name])),
                                  transcription_factor=module.transcription_factor,
                                  gene2weights=leading_edge(rcc, avg2stdrcc, ranking, genes, module.noweights())))

    # Aggregate these regulomes into a single one using the union operator.
    return reduce(Regulome.union, regulomes)


# The 'bincount' implementation is slower but can generated weighted recovery curves.
module2regulome = module2regulome_bincount_impl


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
    # Credits to Thomas Moerman (arboretum package).

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
    def __init__(self, db: Type[RankingDatabase], modules: Sequence[Regulome],
                 motif_annotations_fname: str, sender,
                 rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                 motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0):
        super().__init__(name=db.name)
        self.database = db
        self.motif_annotations_fname = motif_annotations_fname
        self.modules = modules
        self.rank_threshold = rank_threshold
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologuous_identity_threshold = orthologuous_identity_threshold
        self.sender = sender


    def run(self):
        # Load ranking database in memory.
        rnkdb = MemoryDecorator(self.database)
        print("Worker for {}: database loaded in memory.".format(self.name))

        # Load motif annotations in memory.
        motif_annotations = load_motif_annotations(self.motif_annotations_fname,
                                                   motif_similarity_fdr=self.motif_similarity_fdr,
                                                   orthologuous_identity_threshold=self.orthologuous_identity_threshold)
        print("Worker for {}: motif annotations loaded in memory.".format(self.name))

        # Apply module2regulome on all modules.
        def module2regulome(module):
            return module2regulome_numba_impl(rnkdb, module, motif_annotations=motif_annotations,
                                              rank_threshold=self.rank_threshold, auc_threshold=self.auc_threshold,
                                              nes_threshold=self.nes_threshold, avgrcc_sample_frac=None)
        is_not_none = lambda r: r is not None
        regulomes = list(filter(is_not_none, map(module2regulome, self.modules)))
        print("Worker for {}: {} regulomes created.".format(self.name, len(regulomes)))

        # Sending information back to parent process.
        self.sender.send(regulomes)
        self.sender.close()
        print("Worker for {}: Done.".format(self.name))


def derive_regulomes(rnkdbs: Sequence[Type[RankingDatabase]], modules: Sequence[Regulome],
                         motif_annotations_fname: str,
                         rank_threshold: int = 1500, auc_threshold: float = 0.05, nes_threshold=3.0,
                         motif_similarity_fdr: float = 0.001, orthologuous_identity_threshold: float = 0.0,
                         avgrcc_sample_frac: float = None,
                         weighted_recovery=False, client_or_address=None,
                         num_workers=None) -> Sequence[Regulome]:
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
    :param num_workers: If not using a cluster, the number of workers to use for the calculation. None of all available CPUs need to be used.
    :param client_or_address: The client of IP address of the scheduler when working with dask.
    :return: A sequence of regulomes.
    """

    # Load motif annotations.
    motif_annotations = load_motif_annotations(motif_annotations_fname,
                                               motif_similarity_fdr=motif_similarity_fdr,
                                               orthologuous_identity_threshold=orthologuous_identity_threshold)

    is_not_none = lambda r: r is not None

    if not client_or_address:
        # This implementation overcomes the I/O-bounded performance by the dask-based parallelized/distributed version.
        # Each worker (subprocess) loads a dedicated ranking database and motif annotation table into its own memory
        # space before consuming module. The implementation of each worker uses the AUC-first numba JIT based implementation
        # of the algorithm.
        assert len(rnkdbs) <= num_workers if num_workers else cpu_count(), "The number of databases is larger than the number of cores."
        print("Using {} workers.".format(len(rnkdbs)))
        receivers = []
        for db in rnkdbs:
            sender, receiver = Pipe()
            receivers.append(receiver)
            Worker(db, modules, motif_annotations_fname, sender).start()
        return reduce(concat, (recv.recv() for recv in receivers))
    else:
        # Create dask graph.
        from cytoolz.curried import filter as filtercur
        dask_graph = delayed(compose(list, filtercur(is_not_none)))(
            (delayed(module2regulome)
                (db, gs, motif_annotations, rank_threshold, auc_threshold, nes_threshold, avgrcc_sample_frac, weighted_recovery)
                    for db in rnkdbs for gs in modules))

        if client_or_address == "local":
            return dask_graph.compute(get=get, num_workers=num_workers if num_workers else cpu_count())
        else:
            # Run via dask.distributed framework.
            #TODO: Problem when using a cluster: Workers are being killed for an unknown reason.
            client, shutdown_callback = _prepare_client(client_or_address)
            try:
                return client.compute(dask_graph)
            finally:
                shutdown_callback(False)

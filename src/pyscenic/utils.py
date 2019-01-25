# -*- coding: utf-8 -*-

import pandas as pd
from urllib.parse import urljoin
from .genesig import Regulon, GeneSignature
from .math import masked_rho4pairs
from itertools import chain
import numpy as np
from functools import partial
from typing import Sequence, Type
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
import logging


LOGGER = logging.getLogger(__name__)


COLUMN_NAME_TF = "TF"
COLUMN_NAME_MOTIF_ID = "MotifID"
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
COLUMN_NAME_REGULATION = "regulation"
COLUMN_NAME_CORRELATION = "rho"
RHO_THRESHOLD = 0.03


def _create_idx_pairs(adjacencies: pd.DataFrame, exp_mtx: pd.DataFrame) -> np.ndarray:
    """
    :precondition: The column index of the exp_mtx should be sorted in ascending order.
            `exp_mtx = exp_mtx.sort_index(axis=1)`
    """

    # Create sorted list of genes that take part in a TF-target link.
    genes = set(adjacencies.TF).union(set(adjacencies.target))
    sorted_genes = sorted(genes)

    # Find column idx in the expression matrix of each gene that takes part in a link. Having the column index of genes
    # sorted as well as the list of link genes makes sure that we can map indexes back to genes! This only works if
    # all genes we are looking for are part of the expression matrix.
    assert len(set(exp_mtx.columns).intersection(genes)) == len(genes)
    symbol2idx = dict(zip(sorted_genes, np.nonzero(exp_mtx.columns.isin(sorted_genes))[0]))

    # Create numpy array of idx pairs.
    return np.array([[symbol2idx[s1], symbol2idx[s2]] for s1, s2 in zip(adjacencies.TF, adjacencies.target)])


def add_correlation(adjacencies: pd.DataFrame, ex_mtx: pd.DataFrame,
                    rho_threshold=RHO_THRESHOLD, mask_dropouts=False) -> pd.DataFrame:
    """
    Add correlation in expression levels between target and factor.

    :param adjacencies: The dataframe with the TF-target links.
    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param rho_threshold: The threshold on the correlation to decide if a target gene is activated
        (rho > `rho_threshold`) or repressed (rho < -`rho_threshold`).
    :param mask_dropouts: Do not use cells in which either the expression of the TF or the target gene is 0 when
        calculating the correlation between a TF-target pair.
    :return: The adjacencies dataframe with an extra column.
    """
    assert rho_threshold > 0, "rho_threshold should be greater than 0."

    # TODO: Use Spearman correlation instead of Pearson correlation coefficient: Using a non-parametric test like
    # Spearman rank correlation makes much more sense because we want to capture monotonic and not specifically linear
    # relationships between TF and target genes.

    # Assessment of best optimization strategy for calculating dropout masked correlations between TF-target expression:
    #
    # Measurement of time performance of masked_rho (with numba JIT): 136 µs ± 932 ns for a single pair of vectors.
    # For a typical dataset this translates into (for a single core):
    # 1. Calculating the rectangular (TFxtarget) correlation matrix:
    #    (1,564 TFs * 19,812 targets * 136 microseconds * 10e-6)/3600.0 ~ 12 hours.
    #    This approach calculates far too much be has the potential for easy parallelization via numba (cf. current
    #    implementation of masked_rho_2d).
    # 2. Calculating only needed TF-target pairs:
    #    (6,732,441 TF-target links * 136 microseconds * 10e-6)/3600.0 ~ 2h 30 mins.
    #    - Many of these gene-gene links will be duplicate so there might be a potential for memoization. However because
    #    the calculation is already quite fast and the memoization would need to take into account the commutativity of
    #    the operation and involves hashing large numerical vectors, the benefit if this memoization might be minimal.
    #    - Calculation of unique pairs already takes substantial amount of time and does not introduce a substantial
    #    reduction in the number of gene-gene pairs to calculate the correlation for: 6,732,441 => 6,630,720 (2 min 9 s).
    #    This is exactly the additional needed for calculating the rho values for these pairs. No gain here.
    #
    # The other options would have been to used the masked array abstraction provided by numpy but this again
    # this not allow for easy parallelization. In addition the corrcoef operation is far slower than the numba
    # JIT implementation: 2.36 ms ± 62 µs per loop.
    #
    # The best combined approach is to calculate rhos for pairs defined by indexes which is the approach implemented
    # below.

    # Calculate Pearson correlation to infer repression or activation.
    if mask_dropouts:
        ex_mtx = ex_mtx.sort_index(axis=1)
        col_idx_pairs = _create_idx_pairs(adjacencies, ex_mtx)
        rhos = masked_rho4pairs(ex_mtx.values, col_idx_pairs, 0.0)
    else:
        genes = list(set(adjacencies[COLUMN_NAME_TF]).union(set(adjacencies[COLUMN_NAME_TARGET])))
        ex_mtx = ex_mtx[ex_mtx.columns[ex_mtx.columns.isin(genes)]]
        corr_mtx = pd.DataFrame(index=ex_mtx.columns, columns=ex_mtx.columns,
                                data=np.corrcoef(ex_mtx.values.T))
        rhos = np.array([corr_mtx[s2][s1] for s1, s2 in zip(adjacencies.TF, adjacencies.target)])

    regulations = (rhos > rho_threshold).astype(int) - (rhos < -rho_threshold).astype(int)
    return pd.DataFrame(data={
        COLUMN_NAME_TF: adjacencies[COLUMN_NAME_TF].values,
        COLUMN_NAME_TARGET: adjacencies[COLUMN_NAME_TARGET].values,
        COLUMN_NAME_WEIGHT: adjacencies[COLUMN_NAME_WEIGHT].values,
        COLUMN_NAME_REGULATION: regulations,
        COLUMN_NAME_CORRELATION: rhos})


def modules4thr(adjacencies, threshold, context=frozenset(), pattern="weight>{:.3f}"):
    """

    :param adjacencies:
    :param threshold:
    :return:
    """
    for tf_name, df_grp in adjacencies[adjacencies[COLUMN_NAME_WEIGHT] > threshold].groupby(by=COLUMN_NAME_TF):
        if len(df_grp) > 0:
            yield Regulon(
                name="Regulon for {}".format(tf_name),
                context=frozenset([pattern.format(threshold)]).union(context),
                transcription_factor=tf_name,
                gene2weight=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


def modules4top_targets(adjacencies, n, context=frozenset()):
    """

    :param adjacencies:
    :param n:
    :return:
    """
    for tf_name, df_grp in adjacencies.groupby(by=COLUMN_NAME_TF):
        module = df_grp.nlargest(n, COLUMN_NAME_WEIGHT)
        if len(module) > 0:
            yield Regulon(
                name="Regulon for {}".format(tf_name),
                context=frozenset(["top{}".format(n)]).union(context),
                transcription_factor=tf_name,
                gene2weight=list(zip(module[COLUMN_NAME_TARGET].values, module[COLUMN_NAME_WEIGHT].values)))


def modules4top_factors(adjacencies, n, context=frozenset()):
    """

    :param adjacencies:
    :param n:
    :return:
    """
    df = adjacencies.groupby(by=COLUMN_NAME_TARGET).apply(lambda grp: grp.nlargest(n, COLUMN_NAME_WEIGHT))
    for tf_name, df_grp in df.groupby(by=COLUMN_NAME_TF):
        if len(df_grp) > 0:
            yield Regulon(
                name=tf_name,
                context=frozenset(["top{}perTarget".format(n)]).union(context),
                transcription_factor=tf_name,
                gene2weight=list(zip(df_grp[COLUMN_NAME_TARGET].values, df_grp[COLUMN_NAME_WEIGHT].values)))


ACTIVATING_MODULE = "activating"
REPRESSING_MODULE = "repressing"


def modules_from_adjacencies(adjacencies: pd.DataFrame,
                             ex_mtx: pd.DataFrame,
                        thresholds=(0.75, 0.90),
                        top_n_targets=(50,),
                        top_n_regulators=(5,10,50),
                        min_genes=20,
                        absolute_thresholds=False,
                        rho_dichotomize=True,
                        keep_only_activating=True,
                        rho_threshold=RHO_THRESHOLD,
                        rho_mask_dropouts=True) -> Sequence[Regulon]:
    """
    Create modules from a dataframe containing weighted adjacencies between a TF and its target genes.

    :param adjacencies: The dataframe with the TF-target links. This dataframe should have the following columns:
        :py:const:`pyscenic.utils.COLUMN_NAME_TF`, :py:const:`pyscenic.utils.COLUMN_NAME_TARGET` and :py:const:`pyscenic.utils.COLUMN_NAME_WEIGHT` .
    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param thresholds: the first method to create the TF-modules based on the best targets for each transcription factor.
    :param top_n_targets: the second method is to select the top targets for a given TF.
    :param top_n_regulators: the alternative way to create the TF-modules is to select the best regulators for each gene.
    :param min_genes: The required minimum number of genes in a resulting module.
    :param absolute_thresholds: Use absolute thresholds or percentiles to define modules based on best targets of a TF.
    :param rho_dichotomize: Differentiate between activating and repressing modules based on the correlation patterns of
        the expression of the TF and its target genes.
    :param keep_only_activating: Keep only modules in which a TF activates its target genes.
    :param rho_threshold: The threshold on the correlation to decide if a target gene is activated
        (rho > `rho_threshold`) or repressed (rho < -`rho_threshold`).
    :param rho_mask_dropouts: Do not use cells in which either the expression of the TF or the target gene is 0 when
        calculating the correlation between a TF-target pair.
    :return: A sequence of regulons.
    """

    # Duplicate genes need to be removed from the expression matrix to avoid lookup problems in the correlation
    # matrix.
    # In addition, also make sure the expression matrix consists of floating point numbers. This requirement might
    # be violated when dealing with raw counts as input.
    ex_mtx = ex_mtx.T[~ex_mtx.columns.duplicated(keep='first')].T.astype(float)

    # To make the pySCENIC code more robust to the selection of the network inference method in the first step of
    # the pipeline, it is better to use percentiles instead of absolute values for the weight thresholds.
    if not absolute_thresholds:
        def iter_modules(adjc, context):
            yield from chain(chain.from_iterable(
                                    modules4thr(adjc, thr, context, pattern="weight>{}%".format(frac*100))
                                    for thr, frac in zip(list(adjacencies[COLUMN_NAME_WEIGHT].quantile(thresholds)), thresholds)),
                             chain.from_iterable(modules4top_targets(adjc, n, context) for n in top_n_targets),
                             chain.from_iterable(modules4top_factors(adjc, n, context) for n in top_n_regulators))
    else:
        def iter_modules(adjc, context):
            yield from chain(chain.from_iterable(modules4thr(adjc, thr, context) for thr in thresholds),
                         chain.from_iterable(modules4top_targets(adjc, n, context) for n in top_n_targets),
                         chain.from_iterable(modules4top_factors(adjc, n, context) for n in top_n_regulators))

    if not rho_dichotomize:
        # Do not differentiate between activating and repressing modules.
        modules_iter = iter_modules(adjacencies, frozenset())
    else:
        # Relationship between TF and its target, i.e. activator or repressor, is derived using the original expression
        # profiles. The Pearson product-moment correlation coefficient is used to derive this information.

        # Add correlation column and create two disjoint set of adjacencies.
        LOGGER.info("Calculating Pearson correlations.")
        adjacencies = add_correlation(adjacencies, ex_mtx,
                                  rho_threshold=rho_threshold, mask_dropouts=rho_mask_dropouts)
        activating_modules = adjacencies[adjacencies[COLUMN_NAME_REGULATION] > 0.0]
        if keep_only_activating:
            modules_iter = iter_modules(activating_modules, frozenset([ACTIVATING_MODULE]))
        else:
            repressing_modules = adjacencies[adjacencies[COLUMN_NAME_REGULATION] < 0.0]
            modules_iter = chain(iter_modules(activating_modules, frozenset([ACTIVATING_MODULE])),
                             iter_modules(repressing_modules, frozenset([REPRESSING_MODULE])))

    # Derive modules for these adjacencies.
    # + Add the transcription factor to the module.
    #   [We are unable to assess if a TF works in a direct self-regulating way, either inhibiting its own expression or
    #    activating it. Therefore the most unbiased way forward is to add the TF to both activating as well as
    #    repressing modules]
    # + Filter for minimum number of genes.
    LOGGER.info("Creating modules.")
    def add_tf(module):
        return module.add(module.transcription_factor)
    return list(filter(lambda m: len(m) >= min_genes, map(add_tf, modules_iter)))


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


COLUMN_NAME_MOTIF_URL = "MotifURL"


def add_motif_url(df: pd.DataFrame, base_url: str):
    """

    :param df:
    :param base_url:
    :return:
    """
    df[("Enrichment", COLUMN_NAME_MOTIF_URL)] = list(map(partial(urljoin, base_url), df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    return df


def load_motifs(fname: str, sep: str = ',') -> pd.DataFrame:
    """

    :param fname:
    :param sep:
    :return:
    """
    from .transform import COLUMN_NAME_TARGET_GENES, COLUMN_NAME_CONTEXT

    df = pd.read_csv(fname, sep=sep, index_col=[0,1], header=[0,1], skipinitialspace=True)
    df[('Enrichment', COLUMN_NAME_CONTEXT)] = df[('Enrichment', COLUMN_NAME_CONTEXT)].apply(lambda s: eval(s))
    df[('Enrichment', COLUMN_NAME_TARGET_GENES)] = df[('Enrichment', COLUMN_NAME_TARGET_GENES)].apply(lambda s: eval(s))
    return df

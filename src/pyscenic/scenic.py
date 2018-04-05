# -*- coding: utf-8 -*-

import argparse
import os
import logging
from dask.diagnostics import ProgressBar
from multiprocessing import cpu_count
from arboretum.algo import grnboost2
from arboretum.utils import load_tf_names

from .utils import load_from_yaml
from .rnkdb import opendb, RankingDatabase
from .prune import prune2df, find_features, _prepare_client
from .aucell import aucell
from .genesig import GeneSignature
from .log import create_logging_handler
import pandas as pd
import sys
import pickle
from typing import Type, Sequence
from .utils import modules_from_adjacencies
from .transform import df2regulons as df2regs


LOGGER = logging.getLogger(__name__)


class NoProgressBar:
    def __enter__(self):
        return self

    def __exit__(*x):
        pass


def _load_modules(fname: str) -> Sequence[Type[GeneSignature]]:
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    if fname.endswith('.yaml') or fname.endswith('.yml'):
        return load_from_yaml(fname)
    elif fname.endswith('.dat'):
        with open(fname, 'rb') as f:
            return pickle.load(f)
    else:
        LOGGER.error("Unknown file format for \"{}\".".format(fname))
        sys.exit(1)


def _load_dbs(fnames: Sequence[str], nomenclature: str) -> Sequence[Type[RankingDatabase]]:
    def get_name(fname):
        return os.path.basename(fname).split(".")[0]
    return [opendb(fname=fname.name, name=get_name(fname.name), nomenclature=nomenclature) for fname in fnames]


FILE_EXTENSION2SEPARATOR = {
    '.tsv': '\t',
    '.csv': ','
}


def _df2modules(args):
    ext = os.path.splitext(args.module_fname.name)[1]
    adjacencies = pd.read_csv(args.module_fname.name, sep=FILE_EXTENSION2SEPARATOR[ext])
    ex_mtx = pd.read_csv(args.expression_mtx_fname, sep='\t', header=0, index_col=0)
    return modules_from_adjacencies(adjacencies, ex_mtx,
                                    nomenclature = args.nomenclature,
                                    thresholds=args.thresholds,
                                    top_n_targets=args.top_n_targets,
                                    top_n_regulators=args.top_n_regulators,
                                    min_genes=args.min_genes)


def _df2regulons(fname, nomenclature):
    ext = os.path.splitext(fname,)[1]
    df = pd.read_csv(fname, sep=FILE_EXTENSION2SEPARATOR[ext], index_col=[0,1], header=[0,1], skipinitialspace=True)
    df[('Enrichment', 'Context')] = df[('Enrichment', 'Context')].apply(lambda s: eval(s))
    df[('Enrichment', 'TargetGenes')] = df[('Enrichment', 'TargetGenes')].apply(lambda s: eval(s))
    return df2regs(df, nomenclature=nomenclature)


def _load_expression_matrix(args):
    ext = os.path.splitext(args.expression_mtx_fname.name)[1]
    if ext not in FILE_EXTENSION2SEPARATOR:
        LOGGER.error("Unknown file format \"{}\"".format(args.expression_mtx_fname.name))
        sys.exit(1)
    ex_mtx = pd.read_csv(args.expression_mtx_fname, sep=FILE_EXTENSION2SEPARATOR[ext], header=0, index_col=0)
    if args.transpose == 'yes':
        ex_mtx = ex_mtx.T
    return ex_mtx


def find_adjacencies_command(args):
    """
    Infer co-expression modules.
    """
    LOGGER.info("Loading expression matrix.")
    ex_mtx = _load_expression_matrix(args)
    tf_names = load_tf_names(args.tfs_fname.name)

    n_total_genes = len(ex_mtx.columns)
    n_matching_genes = len(ex_mtx.columns.isin(tf_names))
    if n_total_genes == 0:
        LOGGER.error("The expression matrix supplied does not contain any genes. Make sure the extension of the file matches the format (tab separation for TSV and comma sepatration for CSV).")
        sys.exit(1)
    if float(n_matching_genes)/n_total_genes < 0.80:
        LOGGER.warning("Expression data is available for less than 80% of the supplied transcription factors.")

    LOGGER.info("Inferring regulatory networks.")
    client, shutdown_callback = _prepare_client(args.client_or_address, num_workers=args.num_workers)
    try:
        network = grnboost2(expression_data=ex_mtx, tf_names=tf_names, verbose=True, client_or_address=client)
    finally:
        shutdown_callback(False)

    LOGGER.info("Writing results to file.")
    network.to_csv(args.output, index=False, sep='\t')


def prune_targets_command(args):
    """
    Prune targets/find enriched features.
    """
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    if any(args.module_fname.name.endswith(ext) for ext in FILE_EXTENSION2SEPARATOR.keys()):
        LOGGER.info("Creating modules.")
        modules = _df2modules(args)
    else:
        LOGGER.info("Loading modules.")
        modules = _load_modules(args.module_fname.name)

    if len(modules) == 0:
        LOGGER.error("Not a single module loaded")
        sys.exit(1)

    LOGGER.info("Loading databases.")
    nomenclature = modules[0].nomenclature
    dbs = _load_dbs(args.database_fname, nomenclature)

    LOGGER.info("Calculating regulons.")
    motif_annotations_fname = args.annotations_fname.name
    calc_func = find_features if args.no_pruning == "yes" else prune2df
    with ProgressBar() if args.mode == "dask_multiprocessing" else NoProgressBar():
        out = calc_func(dbs, modules, motif_annotations_fname,
                           rank_threshold=args.rank_threshold,
                           auc_threshold=args.auc_threshold,
                           nes_threshold=args.nes_threshold,
                           client_or_address=args.mode,
                           module_chunksize=args.chunk_size,
                           num_workers=args.num_workers)

    LOGGER.info("Writing results to file.")
    out.to_csv(args.output)


def aucell_command(args):
    """
    Calculate regulon enrichment (as AUC values) for cells.
    """
    LOGGER.info("Loading expression matrix.")
    ex_mtx = _load_expression_matrix(args)

    if any(args.regulons_fname.name.endswith(ext) for ext in FILE_EXTENSION2SEPARATOR.keys()):
        LOGGER.info("Creating regulons.")
        regulons = _df2regulons(args.regulons_fname.name, args.nomenclature)
    elif args.regulons_fname.name.endswith('.gmt'):
        LOGGER.info("Loading regulons.")
        regulons = GeneSignature.from_gmt(args.regulons_fname.name, args.nomenclature,
                                           field_separator='\t', gene_separator='\t')
    else:
        LOGGER.info("Loading regulons.")
        regulons = _load_modules(args.regulons_fname.name)

    LOGGER.info("Calculating enrichment.")
    auc_heatmap = aucell(ex_mtx, regulons, auc_threshold=args.auc_threshold,
                         noweights=args.weights != 'yes', num_cores=args.num_workers)

    LOGGER.info("Writing results to file.")
    auc_heatmap.to_csv(args.output)


def add_recovery_parameters(parser):
    group = parser.add_argument_group('motif enrichment arguments')
    group.add_argument('--rank_threshold',
                       type=int, default=5000,
                       help='The rank threshold used for deriving the target genes of an enriched motif (default: 5000).')
    group.add_argument('--auc_threshold',
                       type=float, default=0.05,
                       help='The threshold used for calculating the AUC of a feature as fraction of ranked genes (default: 0.05).')
    group.add_argument('--nes_threshold',
                       type=float, default=3.0,
                       help='The Normalized Enrichment Score (NES) threshold for finding enriched features (default: 3.0).')
    return parser


def add_annotation_parameters(parser):
    group = parser.add_argument_group('motif annotation arguments')
    group.add_argument('--min_orthologous_identity',
                       type=float, default=0.0,
                       help='Minimum orthologous identity to use when annotating enriched motifs (default: 0.0).')
    group.add_argument('--max_similarity_fdr',
                       type=float, default=0.001,
                       help='Maximum FDR in motif similarity to use when annotating enriched motifs (default: 0.001).')
    group.add_argument('--annotations_fname',
                       type=argparse.FileType('r'),
                       help='The name of the file that contains the motif annotations to use.', required=True)
    return parser


def add_module_parameters(parser):
    group = parser.add_argument_group('module generation arguments')
    group.add_argument('--thresholds',
                       type=float, nargs='+', default=[0.001,0.005],
                       help='The first method to create the TF-modules based on the best targets for each transcription factor (default: 0.001 0.005).')
    group.add_argument('--top_n_targets',
                       type=int, nargs='+', default=[50],
                       help='The second method is to select the top targets for a given TF. (default: 50)')
    group.add_argument('--top_n_regulators',
                       type=int, nargs='+', default=[5,10,50],
                       help='The alternative way to create the TF-modules is to select the best regulators for each gene. (default: 5 10 50)')
    group.add_argument('--min_genes',
                       type=int, default=20,
                       help='The minimum number of genes in a module (default: 20).')
    group.add_argument('--nomenclature',
                       type=str, default="HGNC",
                       help='The nomenclature (default: HGNC).')
    group.add_argument('--expression_mtx_fname',
                       type=argparse.FileType('r'),
                       help='The name of the file that contains the expression matrix (CSV).'
                            ' (Only required if modules need to be generated)')
    return parser


def add_computation_parameters(parser):
    group = parser.add_argument_group('computation arguments')
    group.add_argument('--num_workers',
                       type=int, default=cpu_count(),
                       help='The number of workers to use. Only valid of using dask_multiprocessing, custom_multiprocessing or local as mode. (default: {}).'.format(cpu_count()))
    group.add_argument('--client_or_address',
                       type=str, default='local',
                       help='The client or the IP address of the dask scheduler to use.'
                            ' (Only required of dask_cluster is selected as mode)')
    return parser


def create_argument_parser():
    parser = argparse.ArgumentParser(prog='pySCENIC',
                                     description='Single-CEll regulatory Network Inference and Clustering',
                                     fromfile_prefix_chars='@', add_help=True,
                                     epilog="Arguments can be read from file using a @args.txt construct.")

    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "grnboost" command
    parser_grn = subparsers.add_parser('grnboost',
                                         help='Derive co-expression modules from expression matrix.')
    parser_grn.add_argument('expression_mtx_fname',
                               type=argparse.FileType('r'),
                               help='The name of the file that contains the expression matrix (CSV; rows=cells x columns=genes).')
    parser_grn.add_argument('tfs_fname',
                               type=argparse.FileType('r'),
                               help='The name of the file that contains the list of transcription factors (TXT; one TF per line).')
    parser_grn.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, i.e. a table of TF-target genes (CSV).')
    parser_grn.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix (rows=genes x columns=cells).')
    add_computation_parameters(parser_grn)
    parser_grn.set_defaults(func=find_adjacencies_command)

    # create the parser for the "ctx" command
    parser_ctx = subparsers.add_parser('ctx',
                                         help='Find enriched motifs for a gene signature and optionally prune targets from this signature based on cis-regulatory cues.')
    parser_ctx.add_argument('module_fname',
                              type=argparse.FileType('r'),
                              help='The name of the file that contains the signature or the co-expression modules (YAML or pickled DAT).'
                                   'A CSV with adjacencies can also be supplied.')
    parser_ctx.add_argument('database_fname',
                              type=argparse.FileType('r'), nargs='+',
                              help='The name(s) of the regulatory feature databases (FEATHER of LEGACY).')
    parser_ctx.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, i.e. a table of enriched motifs and target genes (CSV).')
    parser_ctx.add_argument('-n', '--no_pruning', action='store_const', const = 'yes',
                              help='Do not perform pruning, i.e. find enriched motifs.')
    parser_ctx.add_argument('--chunk_size',
                       type=int, default=1,
                       help='The size of the module chunks assigned to a node in the dask graph (default: 1).')
    parser_ctx.add_argument('--mode',
                       choices=['custom_multiprocessing', 'dask_multiprocessing', 'dask_cluster'],
                       default='dask_multiprocessing',
                       help='The mode to be used for computing (default: dask_multiprocessing).')
    add_recovery_parameters(parser_ctx)
    add_annotation_parameters(parser_ctx)
    add_computation_parameters(parser_ctx)
    add_module_parameters(parser_ctx)
    parser_ctx.set_defaults(func=prune_targets_command)

    # create the parser for the "aucell" command
    parser_aucell = subparsers.add_parser('aucell', help='Find enrichment of regulons across single cells.')
    parser_aucell.add_argument('expression_mtx_fname',
                            type=argparse.FileType('r'),
                            help='The name of the file that contains the expression matrix (CSV; rows=cells x columns=genes).')
    parser_aucell.add_argument('regulons_fname',
                          type=argparse.FileType('r'),
                               help='The name of the file that contains the co-expression modules (YAML or pickled DAT).'
                                    'A CSV with adjacencies can also be supplied or a GMT file containing gene signatures.')
    parser_aucell.add_argument('-n','--nomenclature',
                               type=str, default='HGNC',
                               help='The nomenclature used for the gene signatures.')
    parser_aucell.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, a matrix of AUC values (CSV; rows=cells x columns=regulons).')
    parser_aucell.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix (rows=genes x columns=cells).')
    parser_aucell.add_argument('-w', '--weights', action='store_const', const = 'yes',
                               help='Use weights associated with genes in recovery analysis.')
    parser_aucell.add_argument('--num_workers',
                       type=int, default=cpu_count(),
                       help='The number of workers to use (default: {}).'.format(cpu_count()))
    add_recovery_parameters(parser_aucell)
    parser_aucell.set_defaults(func=aucell_command)

    return parser


def scenic(argv=None):
    # Set logging level.
    logging_debug_opt = False
    LOGGER.addHandler(create_logging_handler(logging_debug_opt))
    LOGGER.setLevel(logging.DEBUG)

    # Parse arguments.
    parser = create_argument_parser()
    args = parser.parse_args(args=argv)
    if not hasattr(args, 'func'):
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    scenic()

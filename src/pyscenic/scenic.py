# -*- coding: utf-8 -*-

import argparse
import os
import logging
from dask.diagnostics import ProgressBar
from multiprocessing import cpu_count
from arboretum.algo import grnboost2
from arboretum.utils import load_tf_names
from .utils import load_from_yaml
from .rnkdb import open as opendb, RankingDatabase
from .prune import prune2df, find_features
from .aucell import aucell
from .genesig import GeneSignature
from .log import create_logging_handler
import pandas as pd
import sys
import pickle
from typing import Type, Sequence
from .utils import modules_from_adjacencies
from .transform import df2regulomes as df2regs


LOGGER = logging.getLogger(__name__)


class NoProgressBar:
    def __enter__(self):
        return self

    def __exit__(*x):
        pass


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
                       help='The number of workers to use (default: {}).'.format(cpu_count()))
    group.add_argument('--chunk_size',
                       type=int, default=100,
                       help='The size of the module chunks assigned to a node in the dask graph (default: 100).')
    group.add_argument('--mode',
                       choices=['custom_multiprocessing', 'dask_multiprocessing', 'dask_cluster'],
                       default='dask_multiprocessing',
                       help='The mode to be used for computing (default: dask_multiprocessing).')
    group.add_argument('--client_or_address',
                       type=str, default='local',
                       help='The client or the IP address of the dask scheduler to use.'
                            ' (Only required of dask_cluster is selected as mode)')
    return parser


def _load_modules(fname: str) -> Sequence[Type[GeneSignature]]:
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    if fname.endswith('.yaml'):
        return load_from_yaml(fname)
    else:
        with open(fname, 'rb') as f:
            return pickle.load(f)


def _load_dbs(fnames: Sequence[str], nomenclature: str) -> Sequence[Type[RankingDatabase]]:
    def get_name(fname):
        return os.path.basename(fname).split(".")[0]
    return [opendb(fname=fname.name, name=get_name(fname.name), nomenclature=nomenclature) for fname in fnames]


def find_adjacencies_command(args):
    LOGGER.info("Loading datasets.")
    ex_matrix = pd.read_csv(args.expression_mtx_fname, sep='\t', header=0, index_col=0)
    if args.transpose == 'yes':
        ex_matrix = ex_matrix.T
    tf_names = load_tf_names(args.tfs_fname.name)

    # Check whether the supplied data
    n_total_genes = len(ex_matrix.columns)
    n_matching_genes = len(ex_matrix.columns.isin(tf_names))
    if float(n_matching_genes)/n_total_genes < 0.80:
        LOGGER.warning("Expression data is available for less than 80% of the supplied transcription factors.")

    LOGGER.info("Calculating co-expression modules.")
    network = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True, client_or_address=args.client_or_address)

    LOGGER.info("Writing results to file.")
    network.to_csv(args.output, index=False, sep='\t')


def find_motifs_command(args):
    LOGGER.info("Loading modules.")
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    if args.module_fname.name.lower().endswith('.gmt'):
        modules = GeneSignature.from_gmt(args.module_fname.name, args.nomenclature)
    else:
        modules = _load_modules(args.module_fname.name)

    LOGGER.info("Loading databases.")
    nomenclature = modules[0].nomenclature
    dbs = _load_dbs(args.database_fname, nomenclature)

    LOGGER.info("Calculating regulomes.")
    motif_annotations_fname = args.annotations_fname.name
    with ProgressBar() if args.mode == "dask_multiprocessing" else NoProgressBar():
        df = find_features(dbs, modules, motif_annotations_fname,
                           rank_threshold=args.rank_threshold,
                           auc_threshold=args.auc_threshold,
                           nes_threshold=args.nes_threshold,
                           client_or_address=args.mode,
                           module_chunksize=args.chunk_size,
                           num_workers=args.num_workers)

    LOGGER.info("Writing results to file.")
    df.to_csv(args.output)


FILE_EXTENSION2SEPARATOR = {
    '.tsv': '\t',
    '.csv': ','
}

def df2modules(args):
    ext = os.path.splitext(args.module_fname.name)[1]
    adjacencies = pd.read_csv(args.module_fname.name, sep=FILE_EXTENSION2SEPARATOR[ext])
    ex_mtx = pd.read_csv(args.expression_mtx_fname, sep='\t', header=0, index_col=0)
    return modules_from_adjacencies(adjacencies, ex_mtx,
                                    nomenclature = args.nomenclature,
                                    thresholds=args.thresholds,
                                    top_n_targets=args.top_n_targets,
                             top_n_regulators=args.top_n_regulators,
                                              min_genes=args.min_genes)


def prune_targets_command(args):
    if any(args.module_fname.name.endswith(ext) for ext in FILE_EXTENSION2SEPARATOR.keys()):
        LOGGER.info("Creating modules.")
        modules = df2modules(args)
    else:
        LOGGER.info("Loading modules.")
        modules = _load_modules(args.module_fname.name)

    LOGGER.info("Loading databases.")
    nomenclature = modules[0].nomenclature
    dbs = _load_dbs(args.database_fname, nomenclature)

    LOGGER.info("Calculating regulomes.")
    motif_annotations_fname = args.annotations_fname.name
    with ProgressBar() if args.mode == "dask_multiprocessing" else NoProgressBar():
        out = prune2df(dbs, modules, motif_annotations_fname,
                           rank_threshold=args.rank_threshold,
                           auc_threshold=args.auc_threshold,
                           nes_threshold=args.nes_threshold,
                           client_or_address=args.mode,
                           module_chunksize=args.chunk_size,
                           num_workers=args.num_workers)

    LOGGER.info("Writing results to file.")
    out.to_csv(args.output)


def df2regulomes(fname, nomenclature):
    ext = os.path.splitext(fname,)[1]
    df = pd.read_csv(fname, sep=FILE_EXTENSION2SEPARATOR[ext], index_col=[0,1], header=[0,1], skipinitialspace=True)
    df[('Enrichment', 'Context')] = df[('Enrichment', 'Context')].apply(lambda s: eval(s))
    df[('Enrichment', 'TargetGenes')] = df[('Enrichment', 'TargetGenes')].apply(lambda s: eval(s))
    return df2regs(df, nomenclature=nomenclature)


def aucell_command(args):
    LOGGER.info("Loading expression matrix.")
    ex_mtx = pd.read_csv(args.expression_mtx_fname, sep='\t', header=0, index_col=0)
    if args.transpose == 'yes':
        ex_mtx = ex_mtx.T

    if any(args.regulomes_fname.name.endswith(ext) for ext in FILE_EXTENSION2SEPARATOR.keys()):
        LOGGER.info("Creating regulomes.")
        regulomes = df2regulomes(args.regulomes_fname.name, args.nomenclature)
    elif args.regulomes_fname.name.endswith('.gmt'):
        LOGGER.info("Loading regulomes.")
        regulomes = GeneSignature.from_gmt(args.regulomes_fname.name, args.nomenclature, field_separator='\t', gene_separator='\t')
    else:
        LOGGER.info("Loading regulomes.")
        regulomes = load_from_yaml(args.regulomes_fname.name)

    LOGGER.info("Calculating enrichment.")
    auc_heatmap = aucell(ex_mtx, regulomes, args.rank_threshold, args.auc_threshold, args.weights != 'yes', args.num_workers)

    LOGGER.info("Writing results to file.")
    auc_heatmap.to_csv(args.output)


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
                               help='The name of the file that contains the expression matrix (CSV).')
    parser_grn.add_argument('tfs_fname',
                               type=argparse.FileType('r'),
                               help='The name of the file that contains the list of transcription factors (TXT).')
    parser_grn.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream.')
    parser_grn.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix.')
    add_computation_parameters(parser_grn)
    parser_grn.set_defaults(func=find_adjacencies_command)

    # create the parser for the "motifs" command
    #TODO: join with prune and have an option for pruning!
    parser_motifs = subparsers.add_parser('motifs',
                                         help='Find enriched motifs for gene signatures.')
    parser_motifs.add_argument('signatures_fname',
                              type=argparse.FileType('r'),
                              help='The name of the file that contains the gene signatures (GMT or YAML).')
    parser_motifs.add_argument('database_fname',
                              type=argparse.FileType('r'), nargs='+',
                              help='The name(s) of the regulatory feature databases (FEATHER of LEGACY).')
    parser_motifs.add_argument('-n','--nomenclature',
                               type=str, default='HGNC',
                               help='The nomenclature used for the gene signatures.')
    parser_motifs.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream.')
    add_recovery_parameters(parser_motifs)
    add_annotation_parameters(parser_motifs)
    add_computation_parameters(parser_motifs)
    parser_motifs.set_defaults(func=find_motifs_command)

    # create the parser for the "prune" command
    parser_prune = subparsers.add_parser('prune',
                                         help='Prune targets from a co-expression module based on cis-regulatory cues.')
    parser_prune.add_argument('module_fname',
                              type=argparse.FileType('r'),
                              help='The name of the file that contains the co-expression modules (YAML or pickled DAT).'
                                   'A TSV with adjacencies can also be supplied.')
    parser_prune.add_argument('database_fname',
                              type=argparse.FileType('r'), nargs='+',
                              help='The name(s) of the regulatory feature databases (FEATHER of LEGACY).')
    parser_prune.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream.')
    add_recovery_parameters(parser_prune)
    add_annotation_parameters(parser_prune)
    add_computation_parameters(parser_prune)
    add_module_parameters(parser_prune)
    parser_prune.set_defaults(func=prune_targets_command)

    # create the parser for the "aucell" command
    parser_aucell = subparsers.add_parser('aucell', help='Find enrichment of regulomes across single cells.')
    parser_aucell.add_argument('expression_mtx_fname',
                            type=argparse.FileType('r'),
                            help='The name of the file that contains the expression matrix (CSV).')
    parser_aucell.add_argument('regulomes_fname',
                          type=argparse.FileType('r'),
                               help='The name of the file that contains the co-expression modules (YAML or pickled DAT).'
                                    'A TSV with adjacencies can also be supplied or a GMT file containing gene signatures.')
    parser_aucell.add_argument('-n','--nomenclature',
                               type=str, default='HGNC',
                               help='The nomenclature used for the gene signatures.')
    parser_aucell.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream.')
    parser_aucell.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix.')
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

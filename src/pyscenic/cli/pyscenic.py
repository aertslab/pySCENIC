# -*- coding: utf-8 -*-

import os

# Set number of threads to use for OpenBLAS.
os.environ["OPENBLAS_NUM_THREADS"] = "1"

# Set number of threads to use for MKL.
os.environ["MKL_NUM_THREADS"] = "1"

import argparse
import logging
from shutil import copyfile
from dask.diagnostics import ProgressBar
from multiprocessing import cpu_count
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

from pyscenic.utils import modules_from_adjacencies
from pyscenic.rnkdb import opendb, RankingDatabase
from pyscenic.prune import prune2df, find_features, _prepare_client
from pyscenic.aucell import aucell
from pyscenic.log import create_logging_handler
import sys
from typing import Type, Sequence
from .utils import load_exp_matrix, load_signatures, save_matrix, save_enriched_motifs, load_adjacencies, load_modules, append_auc_mtx, ATTRIBUTE_NAME_CELL_IDENTIFIER, ATTRIBUTE_NAME_GENE

try:
    from pyscenic._version import get_versions
    VERSION = get_versions()['version']
except:
    VERSION = "?.?.?"

LOGGER = logging.getLogger(__name__)


def find_adjacencies_command(args):
    """
    Infer co-expression modules.
    """
    LOGGER.info("Loading expression matrix.")
    try:
        ex_mtx = load_exp_matrix(args.expression_mtx_fname.name,
                                 (args.transpose == 'yes'),
                                 args.cell_id_attribute,
                                 args.gene_attribute)
    except ValueError as e:
        LOGGER.error(e)
        sys.exit(1)

    tf_names = load_tf_names(args.tfs_fname.name)

    n_total_genes = len(ex_mtx.columns)
    n_matching_genes = len(ex_mtx.columns.isin(tf_names))
    if n_total_genes == 0:
        LOGGER.error("The expression matrix supplied does not contain any genes. "
                     "Make sure the extension of the file matches the format (tab separation for TSV and "
                     "comma sepatration for CSV).")
        sys.exit(1)
    if float(n_matching_genes)/n_total_genes < 0.80:
        LOGGER.warning("Expression data is available for less than 80% of the supplied transcription factors.")

    LOGGER.info("Inferring regulatory networks.")
    client, shutdown_callback = _prepare_client(args.client_or_address, num_workers=args.num_workers)
    method = grnboost2 if args.method == 'grnboost2' else genie3
    try:
        network = method(expression_data=ex_mtx, tf_names=tf_names, verbose=True, client_or_address=client)
    finally:
        shutdown_callback(False)

    LOGGER.info("Writing results to file.")

    extension = os.path.splitext(args.output.name)[1].lower()
    separator = '\t' if extension == 'tsv' else ','
    network.to_csv(args.output, index=False, sep=separator)


def adjacencies2modules(args):
    try:
        adjacencies = load_adjacencies(args.module_fname.name)
    except ValueError as e:
        LOGGER.error(e)
        sys.exit(1)

    LOGGER.info("Loading expression matrix.")
    try:
        ex_mtx = load_exp_matrix(args.expression_mtx_fname.name,
                                 (args.transpose == 'yes'),
                                 args.cell_id_attribute,
                                 args.gene_attribute)
    except ValueError as e:
        LOGGER.error(e)
        sys.exit(1)

    return modules_from_adjacencies(adjacencies,
                                    ex_mtx,
                                    thresholds=args.thresholds,
                                    top_n_targets=args.top_n_targets,
                                    top_n_regulators=args.top_n_regulators,
                                    min_genes=args.min_genes,
                                    keep_only_activating=(args.all_modules != "yes"))


def _load_dbs(fnames: Sequence[str]) -> Sequence[Type[RankingDatabase]]:
    def get_name(fname):
        return os.path.basename(fname).split(".")[0]
    return [opendb(fname=fname.name, name=get_name(fname.name)) for fname in fnames]


class NoProgressBar:
    def __enter__(self):
        return self

    def __exit__(*x):
        pass


def prune_targets_command(args):
    """
    Prune targets/find enriched features.
    """
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    # The alternative for which was opted in the end is binary pickling.
    extension = os.path.splitext(args.module_fname.name)[1].lower()
    if extension in {'.csv', '.tsv'}:
        if args.expression_mtx_fname is None:
            LOGGER.error("No expression matrix is supplied.")
            sys.exit(0)
        LOGGER.info("Creating modules.")
        modules = adjacencies2modules(args)
    else:
        LOGGER.info("Loading modules.")
        try:
            modules = load_modules(args.module_fname.name)
        except ValueError as e:
            LOGGER.error(e)
            sys.exit(1)

    if len(modules) == 0:
        LOGGER.error("Not a single module loaded")
        sys.exit(1)

    LOGGER.info("Loading databases.")
    dbs = _load_dbs(args.database_fname)

    LOGGER.info("Calculating regulons.")
    motif_annotations_fname = args.annotations_fname.name
    calc_func = find_features if args.no_pruning == "yes" else prune2df
    with ProgressBar() if args.mode == "dask_multiprocessing" else NoProgressBar():
        df_motifs = calc_func(dbs, modules, motif_annotations_fname,
                           rank_threshold=args.rank_threshold,
                           auc_threshold=args.auc_threshold,
                           nes_threshold=args.nes_threshold,
                           client_or_address=args.mode,
                           module_chunksize=args.chunk_size,
                           num_workers=args.num_workers)

    LOGGER.info("Writing results to file.")
    if args.output.name == '<stdout>':
        df_motifs.to_csv(args.output)
    else:
        save_enriched_motifs(df_motifs, args.output.name)


def aucell_command(args):
    """
    Calculate regulon enrichment (as AUC values) for cells.
    """
    LOGGER.info("Loading expression matrix.")
    try:
        ex_mtx = load_exp_matrix(args.expression_mtx_fname.name,
                                 (args.transpose == 'yes'),
                                 args.cell_id_attribute,
                                 args.gene_attribute)
    except ValueError as e:
        LOGGER.error(e)
        sys.exit(1)

    LOGGER.info("Loading gene signatures.")
    try:
        signatures = load_signatures(args.signatures_fname.name)
    except ValueError as e:
        LOGGER.error(e)
        sys.exit(1)

    LOGGER.info("Calculating cellular enrichment.")
    auc_mtx = aucell(ex_mtx, signatures,
                         auc_threshold=args.auc_threshold,
                         noweights=(args.weights != 'yes'),
                         num_workers=args.num_workers)

    LOGGER.info("Writing results to file.")
    extension = os.path.splitext(args.output.name)[1].lower()
    if extension == '.loom':
        try:
            copyfile(args.expression_mtx_fname.name, args.output.name)
            append_auc_mtx(args.output.name, auc_mtx, signatures)
        except OSError as e:
            LOGGER.error("Expression matrix should be provided in the loom file format.")
            sys.exit(1)
    elif args.output.name == '<stdout>':
        transpose = (args.transpose == 'yes')
        (auc_mtx.T if transpose else auc_mtx).to_csv(args.output)
    else:
        save_matrix(auc_mtx, args.output.name, (args.transpose == 'yes'))



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
                       type=float, nargs='+', default=[0.75,0.90],
                       help='The first method to create the TF-modules based on the best targets for each transcription factor (default: 0.75 0.90).')
    group.add_argument('--top_n_targets',
                       type=int, nargs='+', default=[50],
                       help='The second method is to select the top targets for a given TF. (default: 50)')
    group.add_argument('--top_n_regulators',
                       type=int, nargs='+', default=[5,10,50],
                       help='The alternative way to create the TF-modules is to select the best regulators for each gene. (default: 5 10 50)')
    group.add_argument('--min_genes',
                       type=int, default=20,
                       help='The minimum number of genes in a module (default: 20).')
    group.add_argument('--expression_mtx_fname',
                       type=argparse.FileType('r'),
                       help='The name of the file that contains the expression matrix for the single cell experiment.'
                            ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).'
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


def add_loom_parameters(parser):
    group = parser.add_argument_group('loom file arguments')
    group.add_argument('--cell_id_attribute',
                       type=str, default=ATTRIBUTE_NAME_CELL_IDENTIFIER,
                       help='The name of the column attribute that specifies the identifiers of the cells in the loom file.')
    group.add_argument('--gene_attribute',
                       type=str, default=ATTRIBUTE_NAME_GENE,
                       help='The name of the row attribute that specifies the gene symbols in the loom file.')
    return parser



def create_argument_parser():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__).split('.')[0],
                                     description='Single-CEll regulatory Network Inference and Clustering ({})'.format(VERSION),
                                     fromfile_prefix_chars='@', add_help=True,
                                     epilog="Arguments can be read from file using a @args.txt construct. "
                                            "For more information on loom file format see http://loompy.org . "
                                            "For more information on gmt file format see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats .")

    subparsers = parser.add_subparsers(help='sub-command help')

    # --------------------------------------------
    # create the parser for the "grn" command
    # --------------------------------------------
    parser_grn = subparsers.add_parser('grn',
                                         help='Derive co-expression modules from expression matrix.')
    parser_grn.add_argument('expression_mtx_fname',
                               type=argparse.FileType('r'),
                            help='The name of the file that contains the expression matrix for the single cell experiment.'
                                 ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).')
    parser_grn.add_argument('tfs_fname',
                               type=argparse.FileType('r'),
                               help='The name of the file that contains the list of transcription factors (TXT; one TF per line).')
    parser_grn.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, i.e. a table of TF-target genes (CSV).')
    parser_grn.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix (rows=genes x columns=cells).')
    parser_grn.add_argument('-m', '--method', choices=['genie3', 'grnboost2'],
                            default='grnboost2',
                            help='The algorithm for gene regulatory network reconstruction (default: grnboost2).')
    add_computation_parameters(parser_grn)
    add_loom_parameters(parser_grn)
    parser_grn.set_defaults(func=find_adjacencies_command)

    #-----------------------------------------
    # create the parser for the "ctx" command
    #-----------------------------------------
    parser_ctx = subparsers.add_parser('ctx',
                                         help='Find enriched motifs for a gene signature and optionally prune targets from this signature based on cis-regulatory cues.')
    parser_ctx.add_argument('module_fname',
                              type=argparse.FileType('r'),
                              help='The name of the file that contains the signature or the co-expression modules. '
                                   'The following formats are supported: CSV or TSV (adjacencies), YAML, GMT and DAT (modules)')
    parser_ctx.add_argument('database_fname',
                              type=argparse.FileType('r'), nargs='+',
                              help='The name(s) of the regulatory feature databases. '
                                   'Two file formats are supported: feather or db (legacy).')
    parser_ctx.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, i.e. a table of enriched motifs and target genes (csv, tsv)'
                                 ' or collection of regulons (yaml, gmt, dat, json).')
    parser_ctx.add_argument('-n', '--no_pruning', action='store_const', const = 'yes',
                              help='Do not perform pruning, i.e. find enriched motifs.')
    parser_ctx.add_argument('--chunk_size',
                       type=int, default=100,
                       help='The size of the module chunks assigned to a node in the dask graph (default: 100).')
    parser_ctx.add_argument('--mode',
                       choices=['custom_multiprocessing', 'dask_multiprocessing', 'dask_cluster'],
                       default='dask_multiprocessing',
                       help='The mode to be used for computing (default: dask_multiprocessing).')
    parser_ctx.add_argument('-a', '--all_modules', action='store_const', const = 'yes', default='no',
                            help='Included positive and negative regulons in the analysis (default: no, i.e. only positive).')
    parser_ctx.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                            help='Transpose the expression matrix (rows=genes x columns=cells).')
    add_recovery_parameters(parser_ctx)
    add_annotation_parameters(parser_ctx)
    add_computation_parameters(parser_ctx)
    add_module_parameters(parser_ctx)
    add_loom_parameters(parser_ctx)
    parser_ctx.set_defaults(func=prune_targets_command)

    #--------------------------------------------
    # create the parser for the "aucell" command
    # -------------------------------------------
    parser_aucell = subparsers.add_parser('aucell',
                                          help='Quantify activity of gene signatures across single cells.')

    # Mandatory arguments
    parser_aucell.add_argument('expression_mtx_fname',
                            type=argparse.FileType('r'),
                            help='The name of the file that contains the expression matrix for the single cell experiment.'
                                 ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).')
    parser_aucell.add_argument('signatures_fname',
                          type=argparse.FileType('r'),
                               help='The name of the file that contains the gene signatures.'
                                    ' Three file formats are supported: gmt, yaml or dat (pickle).')
    # Optional arguments
    parser_aucell.add_argument('-o', '--output',
                            type=argparse.FileType('w'), default=sys.stdout,
                            help='Output file/stream, a matrix of AUC values.'
                                 ' Two file formats are supported: csv or loom.'
                                 ' If loom file is specified the loom file while contain the original expression matrix and the'
                                 ' calculated AUC values as extra column attributes.')
    parser_aucell.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                               help='Transpose the expression matrix if supplied as csv (rows=genes x columns=cells).')
    parser_aucell.add_argument('-w', '--weights', action='store_const', const = 'yes',
                               help='Use weights associated with genes in recovery analysis.'
                                    ' Is only relevant when gene signatures are supplied as json format.')
    parser_aucell.add_argument('--num_workers',
                       type=int, default=cpu_count(),
                       help='The number of workers to use (default: {}).'.format(cpu_count()))
    add_recovery_parameters(parser_aucell)
    add_loom_parameters(parser_aucell)
    parser_aucell.set_defaults(func=aucell_command)

    return parser


def main(argv=None):
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
    main()

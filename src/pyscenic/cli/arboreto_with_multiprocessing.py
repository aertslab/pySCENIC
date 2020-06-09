#!/usr/bin/env python3

import sys
import time
import loompy as lp
import pandas as pd
from multiprocessing import Pool, cpu_count
import argparse
import tqdm

from arboreto.utils import load_tf_names
from arboreto.algo import genie3, grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network

from pyscenic.cli.utils import load_exp_matrix, suffixes_to_separator
from pathlib import PurePath


################################################################################
################################################################################

parser_grn = argparse.ArgumentParser(description='Run Arboreto using a multiprocessing pool')

parser_grn.add_argument('expression_mtx_fname',
        type=str,
        help='The name of the file that contains the expression matrix for the single cell experiment.'
        ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).')
parser_grn.add_argument('tfs_fname',
        type=str,
        help='The name of the file that contains the list of transcription factors (TXT; one TF per line).')
parser_grn.add_argument('-m', '--method', choices=['genie3', 'grnboost2'],
        default='grnboost2',
        help='The algorithm for gene regulatory network reconstruction (default: grnboost2).')
parser_grn.add_argument('-o', '--output',
        type=str, default=sys.stdout,
        help='Output file/stream, i.e. a table of TF-target genes (TSV).')
parser_grn.add_argument('--num_workers',
        type=int, default=cpu_count(),
        help='The number of workers to use. (default: {}).'.format(cpu_count()))
parser_grn.add_argument('--seed', type=int, required=False, default=None,
        help='Seed value for regressor random state initialization (optional)')

parser_grn.add_argument('--cell_id_attribute',
        type=str, default='CellID',
        help='The name of the column attribute that specifies the identifiers of the cells in the loom file.')
parser_grn.add_argument('--gene_attribute',
        type=str, default='Gene',
        help='The name of the row attribute that specifies the gene symbols in the loom file.')
parser_grn.add_argument('--sparse', action='store_const', const=True, default=False,
        help='If set, load the expression data as a sparse (CSC) matrix.')
parser_grn.add_argument('-t', '--transpose', action='store_const', const = 'yes',
        help='Transpose the expression matrix (rows=genes x columns=cells).')

args = parser_grn.parse_args()


################################################################################
################################################################################


if(args.method == 'grnboost2'):
    method_params = [
        'GBM',      # regressor_type
        SGBM_KWARGS # regressor_kwargs
        ]
elif(args.method == 'genie3'):
    method_params = [
        'RF',       # regressor_type
        RF_KWARGS   # regressor_kwargs
        ]


def run_infer_partial_network(target_gene_index):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=args.seed)
    return( n )


if __name__ == '__main__':

    start_time = time.time()
    ex_matrix = load_exp_matrix(args.expression_mtx_fname,
                             (args.transpose == 'yes'),
                             args.sparse,
                             args.cell_id_attribute,
                             args.gene_attribute)

    if args.sparse:
        gene_names = ex_matrix[1]
        ex_matrix = ex_matrix[0]
    else:
        gene_names = ex_matrix.columns

    end_time = time.time()
    print(f'Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time} seconds...', file=sys.stdout)

    tf_names = load_tf_names(args.tfs_fname)
    print(f'Loaded {len(tf_names)} TFs...', file=sys.stdout)

    ex_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_names)

    print(f'starting {args.method} using {args.num_workers} processes...', file=sys.stdout)
    start_time = time.time()

    with Pool(args.num_workers) as p:
        adjs = list(tqdm.tqdm(p.imap(run_infer_partial_network,
                                     target_gene_indices(gene_names, target_genes='all'),
                                     chunksize=1
                                     ),
                              total=len(gene_names)))

    adj = pd.concat(adjs).sort_values(by='importance', ascending=False)

    end_time = time.time()
    print(f'Done in {end_time - start_time} seconds.', file=sys.stdout)

    extension = PurePath(args.output).suffixes
    adj.to_csv(args.output, index=False, sep=suffixes_to_separator(extension))


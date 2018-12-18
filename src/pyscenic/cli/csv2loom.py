# -*- coding: utf-8 -*-

import os
import argparse
from pyscenic.cli.utils import save_df_as_loom, load_exp_matrix


def create_argument_parser():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__).split('.')[0],
                                     description="Convert a expression matrix from csv to loom file format.",
                                     add_help=True)
    parser.add_argument('csv_fname',
                        type=argparse.FileType('rt'),
                        help='The name of the expression matrix to read in csv format (rows=cells x columns=genes).')
    parser.add_argument('loom_fname',
                        type=argparse.FileType('wb'),
                        help='The name of the expression matrix to create (loom format).')
    parser.add_argument('-t', '--transpose', action='store_const', const = 'yes',
                           help='Transpose the expression matrix (rows=genes x columns=cells).')
    return parser


def convert(fname_csv, fname_loom, transpose):
    df = load_exp_matrix(fname_csv, transpose)
    save_df_as_loom(df, fname_loom)


def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    convert(args.csv_fname.name, args.loom_fname.name, (args.transpose == 'yes'))


if __name__ == "__main__":
    main()

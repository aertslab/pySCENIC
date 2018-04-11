# -*- coding: utf-8 -*-

import os
import argparse
from pyscenic.rnkdb import opendb, InvertedRankingDatabase


def derive_db_name(fname:str) -> str:
    return os.path.basename(fname).split(".")[0]


def create_argument_parser():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__).split('.')[0],
                                     description="Inverts a rankings database to reduce disk volume.",
                                     fromfile_prefix_chars='@', add_help=True)
    parser.add_argument('db_fnames', nargs='+',
                        type=argparse.FileType('rb'),
                        help='The name of the databases in legacy SQL or feather format.')
    parser.add_argument('-o', '--outputdir',
                        type=str, default=os.getcwd(),
                        help='Output directory (default: current directory).')
    parser.add_argument('-n', '--topn',
                        type=int, default=50000,
                        help='The number of top genes/regions to keep in database (default: 50k).')
    return parser


def convert(out_folder, in_fnames, topn):
    for fname in in_fnames:
        print("Inverting {}".format(fname.name))
        name = derive_db_name(fname.name)
        InvertedRankingDatabase.invert(opendb(fname=fname.name, name=name),
                                       os.path.join(out_folder, "{}.inverted.feather".format(name)),
                                       topn)


def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    if len(args.db_fnames) == 0:
        parser.print_help()
    else:
        convert(args.outputdir, args.db_fnames, args.topn)


if __name__ == "__main__":
    main()

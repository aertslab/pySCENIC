# -*- coding: utf-8 -*-

import os
import argparse
from pyscenic.rnkdb import convert2feather


def derive_db_name(fname:str) -> str:
    return os.path.basename(fname).split(".")[0]


def create_argument_parser():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__).split('.')[0],
                                     description="Convert a rankings database in legacy SQL format to the new feather format.",
                                     fromfile_prefix_chars='@', add_help=True)
    parser.add_argument('db_fnames', nargs='+',
                        type=argparse.FileType('rb'),
                        help='The name of the databases in legacy SQL format.')
    parser.add_argument('-o', '--outputdir',
                        type=str, default=os.getcwd(),
                        help='Output directory (default: current directory).')
    return parser


def convert(out_folder, in_fnames):
    for fname in in_fnames:
        print("Converting {}".format(fname.name))
        convert2feather(fname.name, out_folder, derive_db_name(fname.name))


def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    if len(args.db_fnames) == 0:
        parser.print_help()
    else:
        convert(args.outputdir, args.db_fnames)


if __name__ == "__main__":
    main()

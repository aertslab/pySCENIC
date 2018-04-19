# coding=utf-8

import os
import argparse
import sys
from pyscenic.genesig import GeneSignature
from pyscenic.regions import RegionRankingDatabase, Delineation, convert


CODE2DELINEATION = {name.lower(): member for name, member in Delineation.__members__.items()}


def create_argument_parser():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__).split('.')[0],
                                     description="Convert a signature of gene symbols to a signature of region identifiers.",
                                     fromfile_prefix_chars='@', add_help=True)
    parser.add_argument('gmt_fname', nargs=1,
                        type=argparse.FileType('rt'),
                        help='The name of the GMT file that contains the signatures.')
    parser.add_argument('db_fname', nargs=1,
                        type=argparse.FileType('rb'),
                        help='The region-based ranking database.')
    parser.add_argument('-d', 'delineation',
                        choices=CODE2DELINEATION.keys(), default=Delineation.HG19_500BP_UP.name.lower(),
                        help='The candidate regulatory region delineation to use (default: {}).'.format(Delineation.HG19_500BP_UP.name.lower()))
    parser.add_argument('-f', '--fraction',
                        type=float, default=0.8,
                        help='The number of top genes/regions to keep in database (default: 0.8).')
    return parser


def gmt2regions(gmt_fname, db_fname, delineation_code, fraction):
    db = RegionRankingDatabase(fname=db_fname, name=os.path.basename(db_fname))
    signatures = GeneSignature.from_gmt(gmt_fname)
    delineation = CODE2DELINEATION[delineation_code]
    for signature in signatures:
        sys.stdout(signature.name + ',' + ','.join(convert(signature, db, delineation, fraction).genes))


def main():
    parser = create_argument_parser()
    args = parser.parse_args()
    if len(args.db_fnames) == 0:
        parser.print_help()
    else:
        gmt2regions(args.gmt_fname, args.db_fname, args.delineation_code, args.fraction)


if __name__ == "__main__":
    main()

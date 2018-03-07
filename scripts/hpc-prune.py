# -*- coding: utf-8 -*-

import os
import sys
import glob
import logging
import argparse
from configparser import ConfigParser
from pyscenic.utils import load_from_yaml
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.prune import prune2df
from dask.diagnostics import ProgressBar
from cytoolz import mapcat


CONFIG_FILENAME = os.path.join(os.path.dirname(__file__), "hpc-prune.ini")


# Setting the root logger for this entire project.
LOGGER = logging.getLogger(__name__.split(".")[0])


def create_logging_handler(debug: bool) -> logging.Handler:
    """
    Create a handler for logging purposes.

    :param debug: Does debug information need to be logged?
    :return: The logging handler.
    """
    # By default stderr is used as the stream for logging.
    ch = logging.StreamHandler(stream=sys.stderr)
    # Logging level DEBUG is less severe than INFO. Therefore when the logging level is set
    # to DEBUG, information will still be outputted. In addition, errors and warnings are more
    # severe than info and therefore will always be outputted to the log.
    ch.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    return ch


class NoProgressBar:
    def __enter__(self):
        return self

    def __exit__(*x):
        pass


def create_argument_parser():
    parser = argparse.ArgumentParser(prog='hpc-prune')
    parser.add_argument('config_filename',
                        type=str, default=CONFIG_FILENAME,
                        help='Name of configuration file to use.')
    parser.add_argument('-i', '--input',
                        type=str,
                        help='Intput file/stream.')
    parser.add_argument('-o', '--output',
                        type=str,
                        help='Output file/stream.')
    return parser


def run(args):
    # Set logging level.
    logging_debug_opt = False
    LOGGER.addHandler(create_logging_handler(logging_debug_opt))
    LOGGER.setLevel(logging.DEBUG)

    LOGGER.info("Using configuration {}.".format(args.config_filename))
    cfg = ConfigParser()
    cfg.read(args.config_filename)

    LOGGER.info("Loading modules.")
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    modules = load_from_yaml(cfg['data']['modules'] if not args.input else args.input)
    # Filter out modules with to few genes.
    min_genes = int(cfg['parameters']['min_genes'])
    modules = list(filter(lambda m: len(m) >= min_genes, modules))

    LOGGER.info("Loading databases.")
    nomenclature = cfg['parameters']['nomenclature']
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    db_fnames = list(mapcat(glob.glob, cfg['data']['databases'].split(";")))
    dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature=nomenclature) for fname in db_fnames]

    LOGGER.info("Calculating regulomes.")
    motif_annotations_fname = cfg['data']['motif_annotations']
    mode= cfg['parameters']['mode']
    with ProgressBar() if mode == "dask_multiprocessing" else NoProgressBar():
        df = prune2df(dbs, modules, motif_annotations_fname,
                                  rank_threshold=int(cfg['parameters']['rank_threshold']),
                                  auc_threshold=float(cfg['parameters']['auc_threshold']),
                                  nes_threshold=float(cfg['parameters']['nes_threshold']),
                                  client_or_address=mode,
                                  module_chunksize=cfg['parameters']['chunk_size'],
                                  num_workers=int(cfg['parameters']['num_cores']))

    LOGGER.info("Writing results to file.")
    df.to_csv(cfg['parameters']['output'] if not args.output else args.output)


if __name__ == "__main__":
    parser = create_argument_parser()
    run(parser.parse_args())

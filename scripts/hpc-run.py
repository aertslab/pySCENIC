# -*- coding: utf-8 -*-

import os
import sys
import glob
import datetime
from configparser import ConfigParser
from pyscenic.utils import load_from_yaml
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.regulome import derive_regulomes
from dask.diagnostics import ProgressBar
from cytoolz import mapcat


CONFIG_FILENAME = os.path.join(os.path.dirname(__file__), "hpc-run.ini")


def run(cfg_fname):
    print("{} - Using configuration {}.".format(datetime.datetime.now(), cfg_fname))
    cfg = ConfigParser()
    cfg.read(cfg_fname)

    print("{} - Loading modules.".format(datetime.datetime.now()))
    # Loading from YAML is extremely slow. Therefore this is a potential performance improvement.
    # Potential improvements are switching to JSON or to use a CLoader:
    # https://stackoverflow.com/questions/27743711/can-i-speedup-yaml
    modules = load_from_yaml(cfg['data']['modules'])

    print("{} - Loading databases.".format(datetime.datetime.now()))
    nomenclature = cfg['parameters']['nomenclature']
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    db_fnames = list(mapcat(glob.glob, cfg['data']['databases'].split(";")))
    dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature=nomenclature) for fname in db_fnames]

    print("{} - Calculating regulomes.".format(datetime.datetime.now()))
    motif_annotations_fname = cfg['data']['motif_annotations']
    mode= cfg['parameters']['mode']
    if mode == "dask_multiprocessing":
        with ProgressBar():
            df = derive_regulomes(dbs, modules, motif_annotations_fname,
                                  rank_threshold=int(cfg['parameters']['rank_threshold']),
                                  auc_threshold=float(cfg['parameters']['auc_threshold']),
                                  nes_threshold=float(cfg['parameters']['nes_threshold']),
                                  output="df",
                                  client_or_address=mode,
                                  module_chunksize=cfg['parameters']['chunk_size'],
                                  num_workers=int(cfg['parameters']['num_cores']))
    else:
        df = derive_regulomes(dbs, modules, motif_annotations_fname,
                              rank_threshold=int(cfg['parameters']['rank_threshold']),
                              auc_threshold=float(cfg['parameters']['auc_threshold']),
                              nes_threshold=float(cfg['parameters']['nes_threshold']),
                              output="df",
                              client_or_address=mode,
                              module_chunksize=cfg['parameters']['chunk_size'],
                              num_workers=int(cfg['parameters']['num_cores']))

    print("{} - Writing results to file.".format(datetime.datetime.now()))
    df.to_csv(cfg['parameters']['output'])


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) == 2 else CONFIG_FILENAME)

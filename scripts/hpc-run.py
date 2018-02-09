# -*- coding: utf-8 -*-

import os
import glob
import datetime
from configparser import ConfigParser
from pyscenic.utils import load_from_yaml
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.regulome import derive_regulomes
from dask.diagnostics import ProgressBar


CONFIG_FILENAME = os.path.join(os.path.dirname(__file__), "hpc-run.ini")


def run():
    cfg = ConfigParser()
    cfg.read(CONFIG_FILENAME)

    print("{} - Loading modules.".format(datetime.datetime.now()))
    modules = load_from_yaml(cfg['data']['modules'])

    print("{} - Loading databases.".format(datetime.datetime.now()))
    db_fnames = glob.glob(cfg['data']['databases'])
    nomenclature = cfg['parameters']['nomenclature']
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature=nomenclature) for fname in db_fnames]

    print("{} - Calculating regulomes.".format(datetime.datetime.now()))
    motif_annotations_fname = cfg['data']['motif_annotations']
    mode= cfg['parameters']['mode']
    if mode == "dask_multiprocessing":
        with ProgressBar():
            df = derive_regulomes(dbs, modules, motif_annotations_fname, output="df", client_or_address=mode)
    else:
        df = derive_regulomes(dbs, modules, motif_annotations_fname, output="df", client_or_address=mode)

    print("{} - Writing results to file.".format(datetime.datetime.now()))
    df.to_csv(cfg['parameters']['output'])


if __name__ == "__main__":
    run()

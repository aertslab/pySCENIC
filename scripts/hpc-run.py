# -*- coding: utf-8 -*-

import os
import glob
from pyscenic.utils import load_from_yaml
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.regulome import derive_regulomes
from dask.diagnostics import ProgressBar


RESOURCES_FOLDER="/user/leuven/304/vsc30402/data/resources/"
DATABASE_FOLDER = "/user/leuven/304/vsc30402/data/databases/"
FEATHER_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
MODULES_FNAME = os.path.join(RESOURCES_FOLDER, "modules_zeisel_2015.yaml")
RESULTS_FNAME = os.path.join(RESOURCES_FOLDER, "regulomes_zeisel_2015.csv")
NOMENCLATURE = "MGI"


def run():
    modules = load_from_yaml(MODULES_FNAME)
    db_fnames = glob.glob(FEATHER_GLOB)
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature="MGI") for fname in db_fnames]
    with ProgressBar():
        df = derive_regulomes(dbs, modules, MOTIF_ANNOTATIONS_FNAME,
                          output="df", client_or_address="dask_multiprocessing")
        df.to_csv(RESULTS_FNAME)


if __name__ == "__main__":
    run()

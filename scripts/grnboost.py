# -*- coding: utf-8 -*-

import pandas as pd
import os, sys, glob
import datetime
from configparser import ConfigParser
from arboretum.algo import grnboost2
from arboretum.utils import load_tf_names
from dask.distributed import LocalCluster, Client
#from cytoolz import mapcat
import logging
import traceback
from functools import partial


# Setting the root logger for this entire project.
LOGGER = logging.getLogger(__name__.split(".")[0])
CONFIG_FILENAME = os.path.join(os.path.dirname(__file__), "grnboost.ini")


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


def process(mtx_fname, tfs, net_fname, client):
    network = grnboost2(expression_data=pd.read_csv(mtx_fname, sep='\t', index_col=0).T, tf_names=tfs, verbose=True, client_or_address=client)
    network.to_csv(net_fname, index=False)


def run(cfg_fname):
    # Read configuration file.
    cfg = ConfigParser()
    cfg.read(cfg_fname)

    # Set logging level.
    logging_debug_opt = cfg["params"]["debug"].lower().strip() in {"yes", "true", "y"}
    LOGGER.addHandler(create_logging_handler(logging_debug_opt))
    LOGGER.setLevel(logging.DEBUG)

    # Derive file names.
    #mtx_fnames = list(mapcat(glob.glob, cfg['data']['mtx_fnames'].split(";")))
    mtx_fnames = glob.glob(cfg['data']['mtx_fnames'])
    tfs = load_tf_names(cfg['data']['tfs_fname'])

    # Derive cluster information.
    not_cluster_ip = 'scheduler_ip' not in cfg['params']
    if not_cluster_ip:
        local_cluster = LocalCluster(n_workers=int(cfg['params']['num_cores']),
                                 threads_per_worker=1)
        client = Client(local_cluster)
    else:
        class DummyClient:
            def close(self):
                pass
        local_cluster = DummyClient()
        client = cfg['params']['scheduler_ip']

    # Remove fnames that already have a corresponding results file.
    def add_output(fname, out_folder):
        basename = os.path.basename(fname).split('.')[0]
        return fname, os.path.join(out_folder, "{}.net.csv".format(basename))
    out_folder = cfg['data']['out_folder']
    for in_fname, out_fname in filter(lambda t: not os.path.exists(t[1]),
                                    map(partial(add_output, out_folder=out_folder),
                                        mtx_fnames)):
        LOGGER.info("Running GRNboost for {}.".format(in_fname))
        try:
            process(in_fname, tfs, out_fname, client)
        except ValueError as e:
            LOGGER.error("Unable to process {} because of \"{}\". Stacktrace:".format(in_fname, str(e)))
            LOGGER.error(traceback.format_exc())

    if not_cluster_ip:
        client.close()
        local_cluster.close()

    print("{} - Done.".format(datetime.datetime.now()))


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) == 2 else CONFIG_FILENAME)

# coding=utf-8

import pandas as pd
import os
import glob
import pickle
from pyscenic.utils import modules_from_adjacencies


RESOURCES_FOLDER="."
NOMENCLATURE = "HGNC"
RHO_THRESHOLDS = [0.03, 0.1, 0.3]
MODULES_EXT = "modules.dat"
EXP_MTX_EXT = "mtx.tsv"
ADJACENCIES_EXT = "net.csv"


def get_name(fname):
    return os.path.basename(fname).split('.')[0]


def exp_mtx_fname(name):
    return os.path.join(RESOURCES_FOLDER, "{}.{}".format(name, EXP_MTX_EXT))


def calc_modules(adjacencies, exp_mtx, name, rho_dichotomize, rho_threshold=None, mask_dropouts=None):
    if rho_dichotomize:
        print('{} - {} masking - rho threshold {}'.format(name, "with" if mask_dropouts else "without", rho_threshold))

        out_fname = os.path.join(RESOURCES_FOLDER,
            "{}.{}.{}".format(name, rho_threshold, MODULES_EXT) if mask_dropouts
            else "{}.{}.no_mask.{}".format(name, rho_threshold, MODULES_EXT))
    else:
        print('{} - all'.format(name))

        out_fname = os.path.join(RESOURCES_FOLDER, "{}.all.{}".format(name, MODULES_EXT))

    if os.path.isfile(out_fname):
        return

    modules = list(modules_from_adjacencies(adjacencies, exp_mtx, NOMENCLATURE,
                                            rho_dichotomize=rho_dichotomize,
                                            rho_threshold=rho_threshold,
                                            rho_mask_dropouts=mask_dropouts))
    print(len(modules))

    with open(out_fname, 'wb') as f:
        pickle.dump(modules, f)


def run():
    for fname in glob.glob(os.path.join(RESOURCES_FOLDER, '*.{}'.format(ADJACENCIES_EXT))):
        name = get_name(fname)
        mtx_fname = exp_mtx_fname(name)

        # Load datasets.
        adjacencies = pd.read_csv(fname)
        exp_mtx = pd.read_csv(mtx_fname, sep='\t', index_col=0).T

        # Calculate modules.
        for rho_threshold in RHO_THRESHOLDS:
            calc_modules(adjacencies, exp_mtx, name, rho_dichotomize=True, rho_threshold=rho_threshold, mask_dropouts=False)
            calc_modules(adjacencies, exp_mtx, name, rho_dichotomize=True, rho_threshold=rho_threshold, mask_dropouts=True)

        calc_modules(adjacencies, exp_mtx, name, rho_dichotomize=False)


if __name__ == "__main__":
    run()

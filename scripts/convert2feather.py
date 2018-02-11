# -*- coding: utf-8 -*-

import os
import sys
from pyscenic.rnkdb import convert2feather


def derive_db_name(fname:str) -> str:
    return os.path.basename(fname).split(".")[0]


def convert(out_folder, nomenclature, *in_fnames):
    for fname in in_fnames:
        print("Converting {}".format(fname))
        convert2feather(fname, out_folder, derive_db_name(fname), nomenclature)


if __name__ == "__main__":
    convert(os.getcwd(), *sys.argv[2:])

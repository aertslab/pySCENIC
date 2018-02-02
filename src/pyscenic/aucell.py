# -*- coding: utf-8 -*-

import pandas as pd
from .recovery import enrichment4cells


def create_rankings(ex_mtx: pd.DataFrame) -> pd.DataFrame:
    """
    Create a whole genome rankings dataframe from a single cell expression profile dataframe.

    :param ex_mtx: The expression profile matrix. The rows should correspond to different cells, the genes to different
        genes.
    :return: A genome rankings dataframe.
    """
    return ex_mtx.rank(axis=1, ascending=False, method='first').astype('int64')


enrichment = enrichment4cells

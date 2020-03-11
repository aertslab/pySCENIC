Frequently Asked Questions
==========================

I am having problems with Dask
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Arboreto package :code:`v0.1.5`, and some steps of the cisTarget step within pySCENIC, seem to depend on an older version of Dask/Distributed.
Using a more recent version of Dask/Distributed can result in some cryptic errors.
It is recommended to use the older version of Dask and Distributed for stability here:

.. code-block:: bash

    pip install dask==1.0.0 distributed'>=1.21.6,<2.0.0'


But in many cases this still results in issues with the GRN step.
An alternative is to use the multiprocessing implementation of Arboreto recently included in pySCENIC (`arboreto_with_multiprocessing.py <https://github.com/aertslab/pySCENIC/blob/master/scripts/arboreto_with_multiprocessing.py>`_).
This script uses the Arboreto and pySCENIC codebase to run GRNBoost2 (or GENIE3) without Dask.
The eliminates the possibility of running the GRN step across multiple nodes, but brings provides additional stability.
The run time is generally equivalent to the Dask implementation using the same number of workers.

The basic usage is:

.. code-block:: bash

    arboreto_with_multiprocessing.py \
        expr_mat.loom \
        allTFs_hg38.txt \
        --method grnboost2 \
        --output adj.tsv \
        --num_workers 20 \
        --seed 777

The algorithm can be selected using the "code:`--method` option (:code:`genie3` or :code:`grnboost2`).
Possible input formats for the expression data are the same as for the pySCENIC CLI: loom, and csv.

Can I create my own ranking databases?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes you can. The code snippet below shows you how to create your own databases:

.. code-block:: python

    from pyscenic.rnkdb import DataFrameRankingDatabase as RankingDatabase
    import numpy as np
    import pandas as pd

    # Every model in a database is represented by a whole genome ranking. The rankings of the genes must be 0-based.
    df = pd.DataFrame(
            data=[[0, 1],
                  [1, 0]],
            index=['Model1', 'Model2'],
            columns=['Symbol1', 'Symbol2'],
            dtype=np.int32)
    RankingDatabase(df, 'custom').save('custom.db')


Can I draw the distribution of AUC values for a regulon across cells?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import pandas as pd
    import matplotlib.pyplot as plt


    def plot_binarization(auc_mtx: pd.DataFrame, regulon_name: str, threshold: float, bins: int=200, ax=None) -> None:
        """
        Plot the "binarization" process for the given regulon.

        :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
        :param regulon_name: The name of the regulon.
        :param bins: The number of bins to use in the AUC histogram.
        :param threshold: The threshold to use for binarization.
        """
        if ax is None:
            ax=plt.gca()
        auc_mtx[regulon_name].hist(bins=bins,ax=ax)

        ylim = ax.get_ylim()
        ax.plot([threshold]*2, ylim, 'r:')
        ax.set_ylim(ylim)
        ax.set_xlabel('AUC')
        ax.set_ylabel('#')
        ax.set_title(regulon_name)


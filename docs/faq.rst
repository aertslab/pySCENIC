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
An alternative is to use the multiprocessing implementation of Arboreto recently included in pySCENIC (`arboreto_with_multiprocessing.py <https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/cli/arboreto_with_multiprocessing.py>`_).
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

The algorithm can be selected using the ``--method`` option (``genie3`` or ``grnboost2``).
Possible input formats for the expression data are the same as for the pySCENIC CLI: loom, and csv.


How can I prioritize the target genes within a particular regulon?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can be useful to have a better idea of which particular target genes within a regulon are more important, especially in the cases of regulons with many genes.
There are multiple possibilities to address this.

1. **iRegulon analysis.** One possibility is to take the pre-refinement modules for the regulon of interest, and export them for analysis in `iRegulon <http://iregulon.aertslab.org/>`_.
   There are six unrefined modules created from the GRN output for each TF of interest, which are generated using multiple approaches.
   The use of iRegulon for the pruning/refinement allows for more control over the pruning of these modules and possibly a better idea of which genes are more important.
   See the last section, *Further exploration of modules directly from the network inference output*, in 
   `this notebook <http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_downstream-analysis.html>`_
   for an example of how to get started.
   A full tutorial on module refinement with iRegulon can be found `here <http://iregulon.aertslab.org/tutorial.html>`_.

2. **pySCENIC multi-runs.** Another approach is to run the whole pySCENIC procedure multiple times (~10-100x).
   Because of the stochastic nature of the GRN step in particular, a slightly varying result is produced for each run, both in terms of the overall regulons found, as well as the target genes for a particular regulon.
   The target genes (and regulons themselves) can then be scored by the number of times it occurs across all runs, and considered as 'high confidence' if they occur in >80% of runs, for example.
   This has the potential to be very computationally intensive for a whole dataset, but if only a few regulons are of interest, pySCENIC can be run using only these (by limiting the TFs included in the TFs file that goes into the GRN step). 
   The multi-runs capability is implemented in the SCENIC section of our `single cell Nextflow pipeline <https://github.com/vib-singlecell-nf/vsn-pipelines>`_.
   The entrypoint `scenic_multiruns <https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#scenic-multiruns-scenic-multiruns-single-sample-scenic-multiruns>`_ provides an automated way to run this procedure.


How can I create a SCope-compatible loom file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pySCENIC is capable of reading and writing loom files.
However, the loom file created in the AUCell step has the pySCENIC results embedded, but no dimensionality reductions are included.
In order to create a loom file that is ready to be uploaded to the `SCope viewer <http://scope.aertslab.org/>`_, we can use a helper script from `VSN Pipelines <https://github.com/vib-singlecell-nf/vsn-pipelines>`_.
These need to be downloaded:

.. code-block:: bash

    wget https://raw.githubusercontent.com/vib-singlecell-nf/scenic/master/bin/add_visualization.py
    wget https://raw.githubusercontent.com/vib-singlecell-nf/scenic/master/bin/export_to_loom.py

The ``add_visualization.py`` script will take as input the loom file created by ``pyscenic aucell`` and add a basic UMAP and t-SNE based on the SCENIC AUCell matrix.
Some additional packages are required for this, in particular ``MulticoreTSNE`` and ``umap``.
The usage is as follows:

.. code-block:: python

    python add_visualization.py \
        --loom_input auc_mtx.loom \
        --loom_output scenic_visualize.loom \
        --num_workers 8

The output loom file is then ready for use in SCope.

This can also be easily run using the Docker image, which already contains all of the necessary packages (it's still necessary to download the above python files, however):

.. code-block:: bash

    docker run -it --rm -v $PWD:$PWD -w $PWD aertslab/pyscenic:0.10.3 \
        python add_visualization.py \
            --loom_input auc_mtx.loom \
            --loom_output scenic_visualize.loom \
            --num_workers 8


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


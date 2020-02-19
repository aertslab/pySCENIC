Tutorial
========

Other
-----

* See pySCENIC notebooks_ with example analyses

Case studies 
------------

These case studies are associated with the SCENICprotocol_ and the resulting loom files can be viewed in the SCope viewer `here <http://scope.aertslab.org/#/Protocol_Cases/Protocol_Cases/>`_.

PBMC 10k dataset (10x Genomics)
_______________________________

Full SCENIC analysis, plus filtering, clustering, visualization, and SCope-ready loom file creation:

  * `Jupyter notebook <https://github.com/aertslab/SCENICprotocol/tree/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb>`_ | `HTML render <http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.html>`_

Extended analysis post-SCENIC:

  * `Jupyter notebook <https://github.com/aertslab/SCENICprotocol/tree/master/notebooks/PBMC10k_downstream-analysis.ipynb>`_ | `HTML render <http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_downstream-analysis.html>`_

Cancer data sets
________________

  * `Jupyter notebook <https://github.com/aertslab/SCENICprotocol/tree/master/notebooks/SCENIC%20Protocol%20-%20Case%20study%20-%20Cancer%20data%20sets.ipynb>`_ | `HTML render <http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/SCENIC%20Protocol%20-%20Case%20study%20-%20Cancer%20data%20sets.html>`_


Zeisel et al. dataset
---------------------

For this tutorial 3,005 single cell transcriptomes taken from the mouse brain (somatosensory cortex and
hippocampal regions) are used as an example [4]_. The analysis is done in a Jupyter_ notebook.

This dataset is also included in the case studies of the SCENIC protocol, and the analysis notebook can be found therein (
`Jupyter notebook <https://github.com/aertslab/SCENICprotocol/tree/master/notebooks/SCENIC%20Protocol%20-%20Case%20study%20-%20Mouse%20brain%20data%20set.ipynb>`_ | 
`HTML render <http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/SCENIC%20Protocol%20-%20Case%20study%20-%20Mouse%20brain%20data%20set.html>`_).

.. caution::
    If you run this from a python script instead of a Jupyter_ notebook, please enclose the code in
    a :code:`if __name__ == '__main__':` construct.


First we import the necessary modules and declare some constants:

.. code-block:: python

    import os
    import glob
    import pickle
    import pandas as pd
    import numpy as np

    from dask.diagnostics import ProgressBar

    from arboreto.utils import load_tf_names
    from arboreto.algo import grnboost2

    from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies, load_motifs
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell

    import seaborn as sns

    DATA_FOLDER="~/tmp"
    RESOURCES_FOLDER="~/resources"
    DATABASE_FOLDER = "~/databases/"
    SCHEDULER="123.122.8.24:8786"
    DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.mc9nr.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')
    SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "GSE60361_C1-3005-Expression.txt")
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")


Preliminary work
~~~~~~~~~~~~~~~~

The scRNA-Seq data is downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361 and loaded into memory:

.. code-block:: python

    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
    ex_matrix.shape

::

    (3005, 19970)

and the list of Transcription Factors (TF) for *Mus musculus* are read from file.
The list of known TFs for Mm was prepared from TFCat (cf. notebooks_ section).

.. code-block:: python

    tf_names = load_tf_names(MM_TFS_FNAME)


Finally the ranking databases are loaded:

.. code-block:: python

    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    dbs

::

        [FeatherRankingDatabase(name="mm9-tss-centered-10kb-10species.mc9nr"),
         FeatherRankingDatabase(name="mm9-500bp-upstream-7species.mc9nr"),
         FeatherRankingDatabase(name="mm9-500bp-upstream-10species.mc9nr"),
         FeatherRankingDatabase(name="mm9-tss-centered-5kb-10species.mc9nr"),
         FeatherRankingDatabase(name="mm9-tss-centered-10kb-7species.mc9nr"),
         FeatherRankingDatabase(name="mm9-tss-centered-5kb-7species.mc9nr")]

Phase I: Inference of co-expression modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the initial phase of the pySCENIC pipeline the single cell expression profiles are used to infer 
co-expression modules from.

Run GENIE3 or GRNBoost2 from arboreto_ to infer co-expression modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The arboreto package is used for this phase of the pipeline. For this notebook only a sample of 1,000 cells is used
for the co-expression module inference is used.

.. code-block:: python

    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

Derive potential regulons from these co-expression modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regulons are derived from adjacencies based on three methods.

The first method to create the TF-modules is to select the best targets for each transcription factor:

1. Targets with importance > the 50th percentile.
2. Targets with importance > the 75th percentile
3. Targets with importance > the 90th percentile.

The second method is to select the top targets for a given TF:

1. Top 50 targets (targets with highest weight)

The alternative way to create the TF-modules is to select the best regulators for each gene (this is actually how GENIE3 internally works). Then, these targets can be assigned back to each TF to form the TF-modules. In this way we will create three more gene-sets:

1. Targets for which the TF is within its top 5 regulators
2. Targets for which the TF is within its top 10 regulators
3. Targets for which the TF is within its top 50 regulators

A distinction is made between modules which contain targets that are being activated and genes that are being repressed. Relationship between TF and its target, i.e. activator or repressor, is derived using the original expression profiles. The Pearson product-moment correlation coefficient is used to derive this information.

In addition, the transcription factor is added to the module and modules that have less than 20 genes are removed.

.. code-block:: python

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))


Phase II: Prune modules for targets with cis regulatory footprints (aka RcisTarget)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

Clusters can be leveraged in the following way:

.. code-block:: python

    # The clusters can be leveraged via the dask framework:
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address=SCHEDULER)

.. caution::
    The nodes of the clusters need to have access to a shared network drive on which the ranking databases are stored.

Reloading the enriched motifs and regulons from file should be done as follows:

.. code-block:: python

    df = load_motifs(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "rb") as f:
        regulons = pickle.load(f)

Phase III: Cellular regulon enrichment matrix (aka AUCell)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We characterize the different cells in a single-cell transcriptomics experiment via the enrichment of the previously discovered
regulons. Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) of the genes that define this regulon.

.. code-block:: python

    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    sns.clustermap(auc_mtx, figsize=(8,8))


References
----------

.. [4] Zeisel, A. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142 (2015).
.. _Jupyter: http://jupyter.org
.. _notebooks: https://github.com/aertslab/pySCENIC/tree/master/notebooks
.. _`SCENICprotocol`: https://github.com/aertslab/SCENICprotocol


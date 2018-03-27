
pySCENIC
========

|buildstatus|_ |pypipackage|_

pySCENIC is a lightning-fast python implementation of the SCENIC_ pipeline (Single-CEll regulatory Network Inference and
Clustering) which enables biologists to infer transcription factors, gene regulatory networks and cell types from
single-cell RNA-seq data.

The pioneering work was done in R and results were published in Nature Methods [1]_.

pySCENIC can be run on a single desktop machine but easily scales to multi-core clusters to analyze thousands of cells
in no time. The latter is achieved via the dask_ framework for distributed computing [2]_.

The pipeline has three steps:

1. First transcription factors (TFs) and their target genes, together defining a regulon, are derived using gene inference methods which solely rely on correlations between expression of genes across cells. The arboretum_ package is used for this step.
2. These regulons are refined by pruning targets that do not have an enrichment for a corresponding motif of the TF effectively separating direct from indirect targets based on the presence of cis-regulatory footprints.
3. Finally, the original cells are differentiated and clustered on the activity of these discovered regulons.

.. note::
The most impactfull speed improvement is introduced by the arboretum_ package in step 1. This package provides an alternative to GENIE3 [3]_ called GRNBoost2. This package can be controlled from within pySCENIC.

.. sidebar:: **Quick Start**

    * `Installation`_
    * `Tutorial`_
    * `Command Line Interface`_
    * See notebooks_
    * Report an issue_
    * Releases at PyPI_

Features
--------

All the functionality of the original R implementation is available and in addition:

1. You can leverage multi-core and multi-node clusters using dask_ and its distributed_ scheduler.
2. We implemented a version of the recovery of input genes that takes into account weights associated with these genes.
3. Regulons, i.e. the regulatory network that connects a TF with its target genes, with targets that are repressed are now also derived and used for cell enrichment analysis.

Installation
------------

The lastest stable release of the **package** itself can be installed via :code:`pip install pyscenic`.

.. caution::
pySCENIC needs a python 3.5 or greater interpreter.

You can also install the bleeding edge (i.e. less stable) version of the package directly from the source:

.. code-block:: bash

    git clone https://github.com/aertslab/pySCENIC.git
    cd pySCENIC/
    pip install .

To successfully use this pipeline you also need **auxilliary datasets**:

1. *Databases ranking the whole genome* of your species of interest based on regulatory features (i.e. transcription factors). Ranking databases are typically stored in the feather_ format.

=================================  ==============  ====================== ============================
  Database                          Species         Search space            # of orthologous species
=================================  ==============  ====================== ============================
hg19-500bp-upstream-10species_      Homo sapiens   [TSS+500bp,TSS[          10
hg19-500bp-upstream-7species_       Homo sapiens   [TSS+500bp,TSS[          7
hg19-tss-centered-10kb-10species_   Homo sapiens   TSS+/-10kbp              10
hg19-tss-centered-10kb-7species_    Homo sapiens   TSS+/-10kbp              7
hg19-tss-centered-5kb-10species_    Homo sapiens   TSS+/-5kbp               10
hg19-tss-centered-5kb-7species_     Homo sapiens   TSS+/-5kbp               7

hg38-10kb-up-and-down-tss_          Homo sapiens   [TSS+10kb,TSS-10kb]      9
hg38-500bp-up-100bp-down-tss_       Homo sapiens   [TSS+500bp,TSS-100bp]    9

mm9-500bp-upstream-10species_       Mus musculus   [TSS+500bp,TSS[          10
mm9-500bp-upstream-7species_        Mus musculus   [TSS+500bp,TSS[          7
mm9-tss-centered-10kb-10species_    Mus musculus   TSS+/-10kbp              10
mm9-tss-centered-10kb-7species_     Mus musculus   TSS+/-10kbp              7
mm9-tss-centered-5kb-10species_     Mus musculus   TSS+/-5kbp               10
mm9-tss-centered-5kb-7species_      Mus musculus   TSS+/-5kbp               7

mm10-10kb-up-and-down-tss_          Mus musculus   [TSS+10kb,TSS-10kb]      9
mm10-500bp-up-100bp-down-tss_       Mus musculus   [TSS+500bp,TSS-100bp]    9

=================================  ==============  ====================== ============================

.. _hg19-500bp-upstream-10species: http://pyscenic.aertslab.org/databases/hg19-500bp-upstream-10species.mc9nr.feather
.. _hg19-500bp-upstream-7species: http://pyscenic.aertslab.org/databases/hg19-500bp-upstream-7species.mc9nr.feather
.. _hg19-tss-centered-10kb-10species: http://pyscenic.aertslab.org/databases/hg19-tss-centered-10kb-10species.mc9nr.feather
.. _hg19-tss-centered-10kb-7species: http://pyscenic.aertslab.org/databases/hg19-tss-centered-10kb-7species.mc9nr.feather
.. _hg19-tss-centered-5kb-10species: http://pyscenic.aertslab.org/databases/hg19-tss-centered-5kb-10species.mc9nr.feather
.. _hg19-tss-centered-5kb-7species: http://pyscenic.aertslab.org/databases/hg19-tss-centered-5kb-7species.mc9nr.feather

.. _hg38-10kb-up-and-down-tss: http://pyscenic.aertslab.org/databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
.. _hg38-500bp-up-100bp-down-tss: http://pyscenic.aertslab.org/databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather

.. _mm9-500bp-upstream-10species: http://pyscenic.aertslab.org/databases/mm9-500bp-upstream-10species.mc9nr.feather
.. _mm9-500bp-upstream-7species: http://pyscenic.aertslab.org/databases/mm9-500bp-upstream-7species.mc9nr.feather
.. _mm9-tss-centered-10kb-10species: http://pyscenic.aertslab.org/databases/mm9-tss-centered-10kb-10species.mc9nr.feather
.. _mm9-tss-centered-10kb-7species: http://pyscenic.aertslab.org/databases/mm9-tss-centered-10kb-7species.mc9nr.feather
.. _mm9-tss-centered-5kb-10species: http://pyscenic.aertslab.org/databases/mm9-tss-centered-5kb-10species.mc9nr.feather
.. _mm9-tss-centered-5kb-7species: http://pyscenic.aertslab.org/databases/mm9-tss-centered-5kb-7species.mc9nr.feather

.. _mm10-10kb-up-and-down-tss: http://pyscenic.aertslab.org/databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
.. _mm10-500bp-up-100bp-down-tss: http://pyscenic.aertslab.org/databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather

2. *Motif annotation* database providing the missing link between an enriched motif and the transcription factor that binds this motif. This pipeline needs a TSV text file where every line represents a particular annotation.

===================  ==============
  Annotations            Species
===================  ==============
`HGNC annotations`_    Homo sapiens
`MGI annotations`_     Mus musculus
===================  ==============

.. _`HGNC annotations`: http://pyscenic.aertslab.org/resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
.. _`MGI annotations`: http://pyscenic.aertslab.org/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl

.. caution::
These ranking databases are 1.1 Gb each so downloading them might take a while. An annotations file is typically 100Mb in size.

Tutorial
--------

For this tutorial 3,005 single cell transcriptomes taken from the mouse brain (somatosensory cortex and
hippocampal regions) are used as an example [4]_. The analysis is done in a Jupyter_ notebook.

First we import the necessary modules and declare some constants:

.. code-block:: python

    import os
    import pandas as pd
    import numpy as np

    from arboretum.utils import load_tf_names
    from arboretum.algo import grnboost2

    from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies, save_to_yaml
    from pyscenic.prune import prune, prune2df
    from pyscenic.aucell import aucell

    import seaborn as sns

    DATA_FOLDER="~/tmp"
    RESOURCES_FOLDER="~/resources"
    DATABASE_FOLDER = "~/databases/"
    FEATHER_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')
    SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "GSE60361_C1-3005-Expression.txt")
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.yaml")
    NOMENCLATURE = "MGI"


Preliminary work
~~~~~~~~~~~~~~~~

The scRNA-Seq data is downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361 and loaded into memory:

.. code-block:: python

    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0)


Subsequently duplicate genes are removed:

.. code-block:: python

    ex_matrix = ex_matrix[~ex_matrix.index.duplicated(keep='first')]
    ex_matrix.shape

::

    (19970, 3005)

and the list of Transcription Factors (TF) for *Mus musculus* are read from file.
The list of known TFs for Mm was prepared from TFCat (cf. notebooks_ section).

.. code-block:: python

    tf_names = load_tf_names(MM_TFS_FNAME)


Finally the ranking databases are loaded:

.. code-block:: python

    db_fnames = glob.glob(FEATHER_GLOB)
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname), nomenclature="MGI") for fname in db_fnames]
    dbs

::

        [FeatherRankingDatabase(name="mm9-tss-centered-10kb-10species",nomenclature=MGI),
         FeatherRankingDatabase(name="mm9-500bp-upstream-7species",nomenclature=MGI),
         FeatherRankingDatabase(name="mm9-500bp-upstream-10species",nomenclature=MGI),
         FeatherRankingDatabase(name="mm9-tss-centered-5kb-10species",nomenclature=MGI),
         FeatherRankingDatabase(name="mm9-tss-centered-10kb-7species",nomenclature=MGI),
         FeatherRankingDatabase(name="mm9-tss-centered-5kb-7species",nomenclature=MGI)]

Phase I: Inference of co-expression modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the initial phase of the pySCENIC pipeline the single cell expression profiles are used to infer
co-expression modules from.

Run GENIE3 or GRNBoost from arboretum_ to infer co-expression modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The arboretum package is used for this phase of the pipeline. For this notebook only a sample of 1,000 cells is used
for the co-expression module inference is used.

.. code-block:: python

    N_SAMPLES = ex_matrix.shape[1] # Full dataset
    adjancencies = grnboost2(expression_data=ex_matrix.T.sample(n=N_SAMPLES, replace=False),
                        tf_names=tf_names, verbose=True)

Derive potential regulons from these co-expression modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regulons are derived from adjacencies based on three methods.

The first method to create the TF-modules is to select the best targets for each transcription factor:

1. Targets with weight > 0.001
2. Targets with weight > 0.005

The second method is to select the top targets for a given TF:

1. Top 50 targets (targets with highest weight)

The alternative way to create the TF-modules is to select the best regulators for each gene (this is actually how GENIE3 internally works). Then, these targets can be assigned back to each TF to form the TF-modules. In this way we will create three more gene-sets:

1. Targets for which the TF is within its top 5 regulators
2. Targets for which the TF is within its top 10 regulators
3. Targets for which the TF is within its top 50 regulators

A distinction is made between modules which contain targets that are being activated and genes that are being repressed. Relationship between TF and its target, i.e. activator or repressor, is derived using the original expression profiles. The Pearson product-moment correlation coefficient is used to derive this information.

In addition, the transcription factor is added to the module and modules that have less than 20 genes are removed.

.. code-block:: python

    modules = list(modules_from_adjacencies(adjacencies, ex_matrix, nomenclature=NOMENCLATURE))


Phase II: Prune modules for targets with cis regulatory footprints (aka RcisTarget)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    regulons = df2regulons(df, NOMENCLATURE)

Directly calculating regulons without the intermediate dataframe of enriched features is also possible:

.. code-block:: python

    regulons = prune(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    save_to_yaml(regulons, REGULONS_FNAME)

Multi-core systems and clusters can leveraged in the following way:

.. code-block:: python

    # The fastest multi-core implementation:
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME,
                        client_or_address="custom_multiprocessing", num_workers=8)
    # or alternatively:
    regulons = prune(dbs, modules, MOTIF_ANNOTATIONS_FNAME,
                        client_or_address="custom_multiprocessing", num_workers=8)

    # The clusters can be leveraged via the dask framework:
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address="local")
    # or alternatively:
    regulons = prune(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address="local")

Phase III: Cellular regulon enrichment matrix (aka AUCell)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We characterize the different cells in a single-cell transcriptomics experiment via the enrichment of the previously discovered
regulons. Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) of the genes that define this regulon.

.. code-block:: python

    auc_mtx = aucell(ex_matrix.T, regulons, num_workers=4)
    sns.clustermap(auc_mtx, figsize=(8,8))

Command Line Interface
----------------------

A command line version of the tool is included. This tool is available after proper installation of the package via :code:`pip`.

.. code-block:: bash

    { ~ }  » pyscenic                                            ~
    usage: pySCENIC [-h] {grnboost,ctx,aucell} ...

    Single-CEll regulatory Network Inference and Clustering

    positional arguments:
      {grnboost,ctx,aucell}
                            sub-command help
        grnboost            Derive co-expression modules from expression matrix.
        ctx                 Find enriched motifs for a gene signature and
                            optionally prune targets from this signature based on
                            cis-regulatory cues.
        aucell              Find enrichment of regulons across single cells.

    optional arguments:
      -h, --help            show this help message and exit

    Arguments can be read from file using a @args.txt construct.

Website
-------

For more information, please visit LCB_ and SCENIC_.

License
-------

GNU General Public License v3

Acknowledgments
---------------

We are grateful to all providers of TF-annotated position weight matrices, in particular Martha Bulyk (UNIPROBE), Wyeth Wasserman and Albin Sandelin (JASPAR), BioBase (TRANSFAC), Scot Wolfe and Michael Brodsky (FlyFactorSurvey) and Timothy Hughes (cisBP).


References
----------

.. [1] Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. Nat Meth 14, 1083–1086 (2017).
.. [2] Rocklin, M. Dask: parallel computation with blocked algorithms and task scheduling. conference.scipy.org
.. [3] Huynh-Thu, V. A. et al. Inferring regulatory networks from expression data using tree-based methods. PLoS ONE 5, (2010).
.. [4] Zeisel, A. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142 (2015).
.. _dask: https://dask.pydata.org/en/latest/
.. _distributed: https://distributed.readthedocs.io/en/latest/
.. _LCB: https://aertslab.org
.. _feather: https://github.com/wesm/feather
.. _arboretum: https://arboretum.readthedocs.io
.. _notebooks: https://github.com/aertslab/pySCENIC/tree/master/notebooks
.. _issue: https://github.com/aertslab/pySCENIC/issues/new
.. _SCENIC: http://scenic.aertslab.org
.. _PyPI: https://pypi.python.org/pypi/pyscenic
.. _Jupyter: http://jupyter.org

.. |buildstatus| image:: https://travis-ci.org/aertslab/pySCENIC.svg?branch=master
.. _buildstatus: https://travis-ci.org/aertslab/pySCENIC

.. |pypipackage| image:: https://badge.fury.io/py/pyscenic.svg
.. _pypipackage: https://badge.fury.io/py/pyscenic

.. |docstatus| image:: https://readthedocs.org/projects/pyscenic/badge/?version=latest
.. _docstatus: http://pyscenic.readthedocs.io/en/latest/?badge=latest


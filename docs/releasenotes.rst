Release Notes
=============

0.12.0 | 2022-08-16
^^^^^^^^^^^^^^^^^^^

* Only databases in Feather v2 format are supported now (`ctxcore <https://github.com/aertslab/ctxcore>`_ ``>= 0.2``),
  which allow uses recent versions of pyarrow (``>=8.0.0``) instead of very old ones (``<0.17``).
  Databases in the new format can be downloaded from https://resources.aertslab.org/cistarget/databases/
  and end with ``*.genes_vs_motifs.rankings.feather`` or ``*.genes_vs_tracks.rankings.feather``.
* Support clustered motif databases.
* Use custom multiprocessing instead of dask, by default.
* Docker image uses python 3.10 and contains only needed pySCENIC dependencies for CLI usage.
* Remove unneeded scripts and notebooks for unused/deprecated database formats.

0.11.2 | 2021-05-07
^^^^^^^^^^^^^^^^^^^

* Split some core cisTarget functions out into a separate repository, `ctxcore <https://github.com/aertslab/ctxcore>`_. This is now a required package for pySCENIC.

0.11.1 | 2021-02-11
^^^^^^^^^^^^^^^^^^^

* Fix bug in motif url construction (#275)
* Fix for export2loom with sparse dataframe (#278)
* Fix sklearn t-SNE import (#285)
* Updates to Docker image (expose port 8787 for Dask dashboard)

0.11.0 | 2021-02-10
^^^^^^^^^^^^^^^^^^^

**Major features:**

* Updated arboreto_ release (GRN inference step) includes:

  * Support for sparse matrices (using the ``--sparse`` flag in ``pyscenic grn``, or passing a sparse matrix to ``grnboost2``/``genie3``).
  * Fixes to avoid dask metadata mismatch error

* Updated cisTarget:

  * Fix for metadata mismatch in ctx prune2df step
  * Support for databases Apache Parquet format
  * Faster loading from feather databases
  * Bugfix: loading genes from a database (previously missing the last gene name in the database)

* Support for Anndata input and output

* Package updates:

  * Upgrade to newer pandas version
  * Upgrade to newer numba version
  * Upgrade to newer versions of dask, distributed

* Input checks and more descriptive error messages.

  * Check that regulons loaded are not empty.

* Bugfixes:

  * In the regulons output from the cisTarget step, the gene weights were incorrectly assigned to their respective target genes (PR #254).
  * Motif url construction fixed when running ctx without pruning
  * Compression of intermediate files in the CLI steps
  * Handle loom files with non-standard gene/cell attribute names
  * Reformat the genesig gmt input/output
  * Fix AUCell output to loom with non-standard loom attributes


0.10.4 | 2020-11-24
^^^^^^^^^^^^^^^^^^^

* Included new CLI option to add correlation information to the GRN adjacencies file. This can be called with ``pyscenic add_cor``.

0.10.3 | 2020-07-15
^^^^^^^^^^^^^^^^^^^

* Integrate arboreto multiprocessing script into pySCENIC CLI
* Skip modules with zero db overlap in cisTarget step
* Additional error message if regulons file is empty
* Additional error if there is a mismatch between the genes present in the GRN and the expression matrix
* Fixed bug in motif url construction when running without pruning


0.10.2 | 2020-06-05
^^^^^^^^^^^^^^^^^^^

* Bugfix for CLI grn step


0.10.1 | 2020-05-17
^^^^^^^^^^^^^^^^^^^

* CLI: file compression (optionally) enabled for intermediate files for the major steps: grn (adjacencies matrix), ctx (regulons), and aucell (auc matrix). Compression is used when the file name argument has a .gz ending.


0.10.3 | 2020-07-15
^^^^^^^^^^^^^^^^^^^

* Integrate arboreto multiprocessing script into pySCENIC CLI
* Skip modules with zero db overlap in cisTarget step
* Additional error message if regulons file is empty
* Additional error if there is a mismatch between the genes present in the GRN and the expression matrix
* Fixed bug in motif url construciton when running without pruning


0.10.2 | 2020-06-05
^^^^^^^^^^^^^^^^^^^

* Bugfix for CLI grn step


0.10.1 | 2020-05-17
^^^^^^^^^^^^^^^^^^^

* CLI: file compression (optionally) enabled for intermediate files for the major steps: grn (adjacencies matrix), ctx (regulons), and aucell (auc matrix). Compression is used when the file name argument has a .gz ending.


0.10.0 | 2020-02-27
^^^^^^^^^^^^^^^^^^^

* Added a helper script `arboreto_with_multiprocessing.py <https://github.com/aertslab/pySCENIC/blob/master/scripts/arboreto_with_multiprocessing.py>`_ that runs the Arboreto GRN algorithms (GRNBoost2, GENIE3) without Dask for compatibility.

* Ability to set a fixed seed in both the AUCell step and in the calculation of regulon thresholds (CLI parameter :code:`--seed`; aucell function parameter :code:`seed`).

* (since 0.9.18) In the modules_from_adjacencies function, the default value of :code:`rho_mask_dropouts` is changed to False. This now matches the behavior of the R version of SCENIC. The cli version has an additional option to turn dropout masking back on (:code:`--mask_dropouts`).



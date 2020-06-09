Release Notes
=============

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



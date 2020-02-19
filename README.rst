pySCENIC
========

|buildstatus|_ |pypipackage|_ |docstatus|_

pySCENIC is a lightning-fast python implementation of the SCENIC_ pipeline (Single-Cell rEgulatory Network Inference and
Clustering) which enables biologists to infer transcription factors, gene regulatory networks and cell types from
single-cell RNA-seq data.

The pioneering work was done in R and results were published in Nature Methods [1]_.

pySCENIC can be run on a single desktop machine but easily scales to multi-core clusters to analyze thousands of cells
in no time. The latter is achieved via the dask_ framework for distributed computing [2]_.

Full documentation available on `Read the Docs <https://pyscenic.readthedocs.io/en/latest/>`_


Overview
--------

The pipeline has three steps:

1. First transcription factors (TFs) and their target genes, together defining a regulon, are derived using gene inference methods which solely rely on correlations between expression of genes across cells. The arboreto_ package is used for this step.
2. These regulons are refined by pruning targets that do not have an enrichment for a corresponding motif of the TF effectively separating direct from indirect targets based on the presence of cis-regulatory footprints.
3. Finally, the original cells are differentiated and clustered on the activity of these discovered regulons.

The most impactful speed improvement is introduced by the arboreto_ package in step 1. This package provides an alternative to GENIE3 [3]_ called GRNBoost2. This package can be controlled from within pySCENIC.


All the functionality of the original R implementation is available and in addition:

1. You can leverage multi-core and multi-node clusters using dask_ and its distributed_ scheduler.
2. We implemented a version of the recovery of input genes that takes into account weights associated with these genes.
3. Regulons, i.e. the regulatory network that connects a TF with its target genes, with targets that are repressed are now also derived and used for cell enrichment analysis.


Website
-------

For more information, please visit LCB_, SCENIC_ (R version), or SCENICprotocol_ (for a Nextflow implementation).


Acknowledgments
---------------

We are grateful to all providers of TF-annotated position weight matrices, in particular Martha Bulyk (UNIPROBE), Wyeth Wasserman and Albin Sandelin (JASPAR), BioBase (TRANSFAC), Scot Wolfe and Michael Brodsky (FlyFactorSurvey) and Timothy Hughes (cisBP).


References
----------

.. [1] Aibar, S. et al. SCENIC: single-cell regulatory network inference and clustering. Nat Meth 14, 1083–1086 (2017).
.. [2] Rocklin, M. Dask: parallel computation with blocked algorithms and task scheduling. conference.scipy.org
.. [3] Huynh-Thu, V. A. et al. Inferring regulatory networks from expression data using tree-based methods. PLoS ONE 5, (2010).
.. [4] Zeisel, A. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142 (2015).

.. |buildstatus| image:: https://travis-ci.org/aertslab/pySCENIC.svg?branch=master
.. _buildstatus: https://travis-ci.org/aertslab/pySCENIC

.. |pypipackage| image:: https://badge.fury.io/py/pyscenic.svg
.. _pypipackage: https://badge.fury.io/py/pyscenic

.. |docstatus| image:: https://readthedocs.org/projects/pyscenic/badge/?version=latest
.. _docstatus: http://pyscenic.readthedocs.io/en/latest/?badge=latest

.. _SCENIC: http://scenic.aertslab.org
.. _dask: https://dask.pydata.org/en/latest/
.. _distributed: https://distributed.readthedocs.io/en/latest/
.. _arboreto: https://arboreto.readthedocs.io
.. _LCB: https://aertslab.org
.. _`SCENICprotocol`: https://github.com/aertslab/SCENICprotocol
.. _notebooks: https://github.com/aertslab/pySCENIC/tree/master/notebooks
.. _issue: https://github.com/aertslab/pySCENIC/issues/new
.. _PyPI: https://pypi.python.org/pypi/pyscenic


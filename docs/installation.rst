Installation and Usage
======================

Installation
------------

Stable
~~~~~~

The latest stable release of the **package** itself can be installed via 

.. code-block:: bash

    pip install pyscenic


.. caution::
    pySCENIC needs a python 3.6 or greater interpreter.


Development
~~~~~~~~~~~

You can also install the bleeding edge (i.e. less stable) version of the package directly from the source:
 
.. code-block:: bash

    git clone https://github.com/aertslab/pySCENIC.git
    cd pySCENIC/
    pip install .


Containers
~~~~~~~~~

**pySCENIC containers** are also available for download and immediate use. In this case, no compiling or installation is required, provided either Docker or Singularity software is installed on the user's system.  Images are available from `Docker Hub`_. Usage of the containers is shown below (`Docker and Singularity Images`_).


Auxiliary datasets
------------------
To successfully use this pipeline you also need **auxilliary datasets**:

1. *Databases ranking the whole genome* of your species of interest based on regulatory features (i.e. transcription factors). Ranking databases are typically stored in the feather_ format and can be downloaded from cisTargetDBs_.
2. *Motif annotation* database providing the missing link between an enriched motif and the transcription factor that binds this motif. This pipeline needs a TSV text file where every line represents a particular annotation.

=======================  ==========================
  Annotations             Species
=======================  ==========================
`HGNC annotations`_       Homo sapiens
`MGI annotations`_        Mus musculus
`Flybase annotations`_    Drosophila melanogaster
=======================  ==========================

.. _`HGNC annotations`: https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
.. _`MGI annotations`: https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
.. _`Flybase annotations`: https://resources.aertslab.org/cistarget/motif2tf/motifs-v8-nr.flybase-m0.001-o0.0.tbl


.. caution::
    These ranking databases are 1.1 Gb each so downloading them might take a while. An annotations file is typically 100Mb in size.

A list of transcription factors is required for the network inference step (GENIE3/GRNBoost2). These lists can be downloaded from `resources section on GitHub <https://github.com/aertslab/pySCENIC/tree/master/resources>`_.


Command Line Interface
----------------------

A command line version of the tool is included. This tool is available after proper installation of the package via :code:`pip`.

.. code-block:: bash

    { ~ }  » pyscenic                                            ~
    usage: pyscenic [-h] {grn,ctx,aucell} ...

    Single-CEll regulatory Network Inference and Clustering (0.9.19)

    positional arguments:
      {grn,ctx,aucell}  sub-command help
        grn             Derive co-expression modules from expression matrix.
        ctx             Find enriched motifs for a gene signature and optionally
                        prune targets from this signature based on cis-regulatory
                        cues.
        aucell          Quantify activity of gene signatures across single cells.

    optional arguments:
      -h, --help        show this help message and exit

    Arguments can be read from file using a @args.txt construct. For more
    information on loom file format see http://loompy.org . For more information
    on gmt file format see https://software.broadinstitute.org/cancer/software/gse
    a/wiki/index.php/Data_formats .


Docker and Singularity Images
-----------------------------

pySCENIC is available to use with both Docker and Singularity, and tool usage from a container is similar to that of the command line interface.
Note that the feather databases, transcription factors, and motif annotation databases need to be accessible to the container via a mounted volume.
In the below examples, a single volume mount is used for simplicity, which will contains the input, output, and databases files.

For additional usage examples, see the documentation associated with the `SCENIC protocol <https://github.com/aertslab/SCENICprotocol/blob/master/docs/installation.md>`_ Nextflow implementation.

Docker
~~~~~~

Docker images are available from `Docker Hub`_, and can be obtained by running :code:`docker pull aertslab/pyscenic:[version]`, with the version tag as the latest release.

To run pySCENIC using Docker, use the following three steps.
A mount point (or more than one) needs to be specified, which contains the input data and necessary resources).

.. code-block:: bash

    docker run -it --rm \
        -v /path/to/data:/scenicdata \
        aertslab/pyscenic:[version] pyscenic grn \
            --num_workers 6 \
            -o /scenicdata/expr_mat.adjacencies.tsv \
            /scenicdata/expr_mat.tsv \
            /scenicdata/allTFs_hg38.txt

    docker run -it --rm \
        -v /path/to/data:/scenicdata \
        aertslab/pyscenic:[version] pyscenic ctx \
            /scenicdata/expr_mat.adjacencies.tsv \
            /scenicdata/hg19-tss-centered-5kb-7species.mc9nr.feather \
            /scenicdata/hg19-tss-centered-10kb-7species.mc9nr.feather \
            --annotations_fname /scenicdata/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
            --expression_mtx_fname /scenicdata/expr_mat.tsv \
            --mode "dask_multiprocessing" \
            --output /scenicdata/regulons.csv \
            --num_workers 6

    docker run -it --rm \
        -v /path/to/data:/scenicdata \
        aertslab/pyscenic:[version] pyscenic aucell \
            /scenicdata/expr_mat.tsv \
            /scenicdata/regulons.csv \
            -o /scenicdata/auc_mtx.csv \
            --num_workers 6

Singularity
~~~~~~~~~~~

As of release :code:`0.9.19`, pySCENIC Singularity images are no longer being built on `Singularity Hub`_, however images can easily be built using Docker Hub as a source:

.. code-block:: bash

    singularity build aertslab-pyscenic-0.9.19.sif docker://aertslab/pyscenic:0.9.19


To run pySCENIC with Singularity, the usage is very similar to that of Docker.
Note that in Singularity 3.0+, the mount points are automatically overlaid, but bind points can be specified similarly to Docker with :code:`--bind`/:code:`-B`.
The first step (GRN inference) is shown as an example:

.. code-block:: bash

    singularity run pySCENIC_0.9.19.sif \
        pyscenic grn \
            --num_workers 6 \
            -o expr_mat.adjacencies.tsv \
            expr_mat.tsv \
            allTFs_hg38.txt


Using the Docker or Singularity images with Jupyter notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As of version 0.9.7, the pySCENIC containers have the `ipykernel` package installed, and can also be used interactively in a notebook.
This can be achieved using a kernel command similar to the following (for singularity).
Note that in this case, a bind needs to be specified.

.. code-block:: bash

    singularity exec -B /data:/data pySCENIC_0.9.7.sif ipython kernel -f {connection_file}


Nextflow
--------

The CLI to pySCENIC has also been streamlined into a pipeline that can be run with a single command, using the Nextflow workflow manager.
For details on this usage, along with more detailed pySCENIC tutorials, see the `SCENICprotocol`_ repository.


.. _`Singularity Hub`: https://www.singularity-hub.org/collections/2033
.. _`SCENICprotocol`: https://github.com/aertslab/SCENICprotocol
.. _dask: https://dask.pydata.org/en/latest/
.. _distributed: https://distributed.readthedocs.io/en/latest/
.. _`Docker Hub`: https://hub.docker.com/r/aertslab/pyscenic
.. _feather: https://github.com/wesm/feather
.. _cisTargetDBs: https://resources.aertslab.org/cistarget/


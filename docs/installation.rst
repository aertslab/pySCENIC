Installation and Usage
======================

Installation
------------

Stable
~~~~~~

The latest stable release of the **package** itself can be installed via 

.. code-block:: bash

    pip install pyscenic


Note that pySCENIC needs some prerequisites installed before running ``pip install`` in a new conda environment.
For example:

.. code-block:: bash

    conda create -y -n pyscenic python=3.10
    conda activate pyscenic

    pip install pyscenic


.. caution::
    pySCENIC needs a python 3.7 or greater interpreter.


Development
~~~~~~~~~~~

You can also install the bleeding edge (i.e. less stable) version of the package directly from the source:
 
.. code-block:: bash

    git clone https://github.com/aertslab/pySCENIC.git
    cd pySCENIC/
    pip install .


Containers
~~~~~~~~~~

**pySCENIC containers** are also available for download and immediate use.
In this case, no compiling or installation is required, provided either Docker/Podman or Singularity/Apptainer software is installed on the user's system.
See (`Docker/Podman and Singularity/Apptainer Images`_).

Auxiliary datasets
------------------

To successfully use this pipeline you also need **auxilliary datasets** available at `cistargetDBs website`_:

1. `Databases ranking the whole genome`_ of your species of interest based on regulatory features (i.e. transcription factors) in `feather`_ format.
2. `Motif to TF annotations`_ database providing the missing link between an enriched motif and the transcription factor that binds this motif. This pipeline needs a TSV text file where every line represents a particular annotation.

.. caution::
    These ranking databases are 1.1 Gb each so downloading them might take a while. An annotations file is typically 100Mb in size.

3. A `list of transcription factors`_ is required for the network inference step (GENIE3/GRNBoost2).


Command Line Interface
----------------------

A command line version of the tool is included. This tool is available after proper installation of the package via :code:`pip`.

.. code-block:: bash

    $ pyscenic -h
    usage: pyscenic [-h] {grn,add_cor,ctx,aucell} ...

    Single-Cell rEgulatory Network Inference and Clustering (0.12.0)

    positional arguments:
      {grn,add_cor,ctx,aucell}
                            sub-command help
        grn                 Derive co-expression modules from expression matrix.
        add_cor             [Optional] Add Pearson correlations based on TF-gene
                            expression to the network adjacencies output from the
                            GRN step, and output these to a new adjacencies file.
                            This will normally be done during the "ctx" step.
        ctx                 Find enriched motifs for a gene signature and
                            optionally prune targets from this signature based on
                            cis-regulatory cues.
        aucell              Quantify activity of gene signatures across single
                            cells.

    optional arguments:
      -h, --help            show this help message and exit

    Arguments can be read from file using a @args.txt construct. For more
    information on loom file format see http://loompy.org . For more information
    on gmt file format see https://software.broadinstitute.org/cancer/software/gse
    a/wiki/index.php/Data_formats .


Docker/Podman and Singularity/Apptainer Images
----------------------------------------------

pySCENIC is available to use with both Docker/Podman and Singularity/Apptainer, and tool usage from a container is similar to that of the command line interface.
Note that the `feather`_ databases, transcription factors, and motif annotation databases need to be accessible to the container via a mounted volume.
In the below examples, a single volume mount is used for simplicity, which will contains the input, output, and databases files.

For additional usage examples, see the documentation associated with the `SCENIC protocol <https://github.com/aertslab/SCENICprotocol/blob/master/docs/installation.md>`_ Nextflow implementation.

Docker/Podman
~~~~~~~~~~~~~

Docker/Podman images are available at `Docker Hub pySCENIC`_ and `Docker Hub pySCENIC with scanpy`_ , and can be obtained by running:

.. code-block:: bash

    # pySCENIC CLI version (recommended).
    docker pull aertslab/pyscenic:0.12.0
    podman pull docker://aertslab/pyscenic:0.12.0

    # pySCENIC CLI version + ipython kernel + scanpy.
    docker pull aertslab/pyscenic_scanpy:0.12.0_1.9.1
    podman pull docker://aertslab/pyscenic_scanpy:0.12.0_1.9.1

To run pySCENIC using Docker/Podman, use the following three steps.
A mount point (or more than one) needs to be specified, which contains the input data and necessary resources).

.. code-block:: bash

    docker run -it --rm \
        -v /data:/data \
        aertslab/pyscenic:0.12.0 pyscenic grn \
            --num_workers 6 \
            -o /data/expr_mat.adjacencies.tsv \
            /data/expr_mat.tsv \
            /data/allTFs_hg38.txt

    docker run -it --rm \
        -v /data:/data \
        aertslab/pyscenic:0.12.0 pyscenic ctx \
            /data/expr_mat.adjacencies.tsv \
            /data/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
            /data/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
            --annotations_fname /data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
            --expression_mtx_fname /data/expr_mat.tsv \
            --mode "custom_multiprocessing" \
            --output /data/regulons.csv \
            --num_workers 6

    docker run -it --rm \
        -v /data:/data \
        aertslab/pyscenic:0.12.0 pyscenic aucell \
            /data/expr_mat.tsv \
            /data/regulons.csv \
            -o /data/auc_mtx.csv \
            --num_workers 6

Singularity/Apptainer
~~~~~~~~~~~~~~~~~~~~~

Singularity/Apptainer images can be build from the Docker Hub image as source:

.. code-block:: bash

    # pySCENIC CLI version.
    singularity build aertslab-pyscenic-0.12.0.sif docker://aertslab/pyscenic:0.12.0
    apptainer build aertslab-pyscenic-0.12.0.sif docker://aertslab/pyscenic:0.12.0

    # pySCENIC CLI version + ipython kernel + scanpy.
    singularity build aertslab-pyscenic-scanpy-0.12.0-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.0_1.9.1
    apptainer build aertslab-pyscenic-0.12.0-1.9.1.sif docker://aertslab/pyscenic_scanpy:0.12.0_1.9.1


To run pySCENIC with Singularity/Apptainer, the usage is very similar to that of Docker/Podman.

.. code-block:: bash

    singularity run aertslab-pyscenic-0.12.0.sif \
        pyscenic grn \
            -B /data:/data
            --num_workers 6 \
            -o /data/expr_mat.adjacencies.tsv \
            /data/expr_mat.tsv \
            /data/allTFs_hg38.txt


Using the Docker/Podman or Singularity/Apptainer images with Jupyter notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pySCENIC containers with scanpy have the ``ipykernel`` package installed, and can also be used interactively in a notebook.
This can be achieved using a kernel command similar to the following (for singularity).
Note that in this case, a bind needs to be specified.

.. code-block:: bash

    singularity exec -B /data:/data aertslab-pyscenic-scanpy-latest.sif ipython kernel -f {connection_file}

More generally, a local or remote kernel can be set up by using the following examples.
These would go in a kernel file in ``~/.local/share/jupyter/kernels/pyscenic-latest/kernel.json`` (for example).

**Remote singularity kernel:**

.. code-block:: bash

    {
      "argv": [
        "/software/jupyter/bin/python",
        "-m",
        "remote_ikernel",
        "--interface",
        "ssh",
        "--host",
        "r23i27n14",
        "--workdir",
        "~/",
        "--kernel_cmd",
        "singularity",
        "exec",
        "-B",
        "/path/to/mounts",
        "/path/to/aertslab-pyscenic-scanpy-latest.sif",
        "ipython",
        "kernel",
        "-f",
        "{connection_file}"
      ],
      "display_name": "pySCENIC singularity remote",
      "language": "Python"
    }

**Local singularity kernel:**

.. code-block:: bash

    {
        "argv": [
         "singularity",
         "exec",
         "-B",
         "/path/to/mounts",
         "/path/to/aertslab-pyscenic-scanpy-latest.sif",
         "ipython",
         "kernel",
         "-f",
         "{connection_file}"
        ],
        "display_name": "pySCENIC singularity local",
        "language": "python"
    }


Nextflow
--------

There are two Nextflow implementations available:

* `SCENICprotocol`_: A Nextflow DSL1 implementation.
* `VSNPipelines`_: A Nextflow DSL2 implementation.


.. _`SCENICprotocol`: https://github.com/aertslab/SCENICprotocol
.. _`VSNPipelines`: https://github.com/vib-singlecell-nf/vsn-pipelines

.. _`Docker Hub pySCENIC`: https://hub.docker.com/r/aertslab/pyscenic
.. _`Docker Hub pySCENIC with scanpy`: https://hub.docker.com/r/aertslab/pyscenic_scanpy

.. _feather: https://arrow.apache.org/docs/python/feather.html

.. _`cistargetDBs website`: https://resources.aertslab.org/cistarget/
.. _`Databases ranking the whole genome`: https://resources.aertslab.org/cistarget/databases/
.. _`Motif to TF annotations`: https://resources.aertslab.org/cistarget/motif2tf/
.. _`list of transcription factors`: https://resources.aertslab.org/cistarget/tf_lists/

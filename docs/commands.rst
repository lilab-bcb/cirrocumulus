Commands
----------------

Launch
^^^^^^^^^^^^^^^

The `launch` command opens one or more datasets in h5ad_, loom_, Seurat_, TileDB_, zarr_, or `STAR-Fusion`_ formats. Results are stored on disk in JSON format.

.. argparse::
   :ref: cirrocumulus.launch.create_parser
   :prog: cirro launch

Serve
^^^^^^^^^^^^^

The `serve` command starts the cirrocumulus server for use in a shared server environment which can handle concurrent requests from multiple users.
The server can optionally enforce permissions at the dataset level, in order to securely share datasets with collaborators.
Additionally, annotations and sets are shared among all users authorized to view a dataset and are stored in a database.
The server can be deployed on a cloud VM, an on-premise machine, or on Google App Engine. When deployed in App Engine, the server can be configured
to be use Google Cloud Firestore as a database. On AWS, we recommend using DocumentDB. Please note that no datasets are available until you import a dataset into cirrocumulus.

.. argparse::
   :ref: cirrocumulus.serve.create_parser
   :prog: cirro serve


Prepare Data
^^^^^^^^^^^^^^


The `prepare_data` command is used to freeze an h5ad_, loom_, or Seurat_ file in cirrocumulus format. The cirrocumulus format allows
efficient partial dataset retrieval over a network (e.g. Google or S3 bucket) using limited memory. Please note that when converting
Seurat_ files, the `data` slot from the default assay is used.


.. argparse::
   :ref: cirrocumulus.prepare_data.create_parser
   :prog: cirro prepare_data


Concatenate
^^^^^^^^^^^^^^
The `concat` command is used to concatenate multiple datasets into a single dataset. Both spatial and non-spatial embeddings are concatenated in a tiled layout:

    .. image:: images/concat.png


.. argparse::
   :ref: cirrocumulus.concat.create_parser
   :prog: cirro concat


.. _h5ad: https://anndata.readthedocs.io/
.. _loom: https://linnarssonlab.org/loompy/format/
.. _STAR-Fusion: https://github.com/STAR-Fusion/STAR-Fusion/wiki
.. _Seurat: https://satijalab.org/seurat/
.. _TileDB: https://tiledb.com/
.. _zarr: https://zarr.readthedocs.io/

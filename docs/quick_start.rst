Quick Start
-------------

Install the package::

    pip install cirrocumulus

Launch cirrocumulus via the command line::

    cirro launch <path_to_dataset>

- Datasets can be provided in h5ad_ or loom_ formats. Additionally cirrocumulus can read in fusions produced by `STAR-Fusion`_.
- Launch accepts more than one dataset to support cases in which modalities are stored in separate files.
- Marker lists can be provided in JSON format to quickly load features from predefined lists.

Run with `example data`_ consisting of 3k PBMCs from a healthy donor::

    cirro launch pbmc3k.h5ad --markers markers.json


Run with `example spatial data`_  from `10x Genomics <https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node/>`_::

    cirro launch V1_Human_Lymph_Node/V1_Human_Lymph_Node.h5ad --spatial V1_Human_Lymph_Node/spatial


.. _example data: https://github.com/klarman-cell-observatory/cirrocumulus/raw/master/docs/example_data.zip
.. _example spatial data: https://github.com/klarman-cell-observatory/cirrocumulus/raw/master/docs/V1_Human_Lymph_Node.zip
.. _h5ad: https://anndata.readthedocs.io/
.. _loom: https://linnarssonlab.org/loompy/format/
.. _STAR-Fusion: https://github.com/STAR-Fusion/STAR-Fusion/wiki

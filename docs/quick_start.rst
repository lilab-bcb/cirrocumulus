Quick Start
-------------

Install the package::

    pip install cirrocumulus

Launch cirrocumulus via the command line::

    cirro launch <path_to_dataset>

- Datasets can be provided in h5ad_, loom_  or `STAR-Fusion`_ format. Seurat objects
  can be loaded after converting to h5ad_ or loom_ format (see vignette_).
- Launch accepts more than one dataset to enable quick dataset switching or to combine modalities (e.g gene fusions and expression) stored in separate files.
- Predefined marker lists can be provided in JSON format to quickly browse features of interest.

Example Data
^^^^^^^^^^^^^
- Download `3k PBMCs from a healthy donor data`_ and launch::

    cirro launch pbmc3k.h5ad --markers markers.json


- Download `human lymph node spatial data`_ and launch::

    cirro launch V1_Human_Lymph_Node/V1_Human_Lymph_Node.h5ad --spatial V1_Human_Lymph_Node/spatial


.. _3k PBMCs from a healthy donor data: https://github.com/klarman-cell-observatory/cirrocumulus/raw/master/docs/example_data.zip
.. _human lymph node spatial data: https://github.com/klarman-cell-observatory/cirrocumulus/raw/master/docs/V1_Human_Lymph_Node.zip
.. _h5ad: https://anndata.readthedocs.io/
.. _loom: https://linnarssonlab.org/loompy/format/
.. _STAR-Fusion: https://github.com/STAR-Fusion/STAR-Fusion/wiki
.. _vignette: https://satijalab.org/seurat/v3.2/conversion_vignette.html

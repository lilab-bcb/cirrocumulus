Quick Start
-------------

Install the package::

    pip install cirrocumulus

Launch cirrocumulus via the command line::

    cirro launch <path_to_dataset>

Run with `example data`_::

    cirro launch pbmc3k.h5ad --markers markers.json

The `example data`_ consists of 3k PBMCs from a Healthy Donor and are freely available from 10x Genomics

- Datasets can be provided in h5ad, loom, or STAR-Fusion formats.
- Launch accepts more than one dataset to support cases in which modalities are stored in separate files.
- Marker lists can be provided in JSON format to quickly load features from predefined lists.


.. _example data: https://github.com/klarman-cell-observatory/cirrocumulus/raw/master/docs/example_data.zip
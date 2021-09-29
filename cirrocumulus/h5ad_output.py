import logging
import os

import pandas._libs.json as ujson

from cirrocumulus.anndata_util import DataType, get_pegasus_marker_keys, get_scanpy_marker_keys

logger = logging.getLogger("cirro")


def save_datasets_h5ad(datasets, schema, output_directory, filesystem, whitelist):
    adata = None
    module_dataset = None
    for dataset in datasets:
        if dataset.uns.get('data_type') == DataType.MODULE:
            module_dataset = dataset
        else:
            adata = dataset
    adata.strings_to_categoricals()
    if module_dataset is not None:
        module_dataset.strings_to_categoricals()
        d = dict(X=module_dataset.X, var=module_dataset.var)
        adata.uns['modules'] = d

    with filesystem.open(os.path.join(output_directory, 'index.json.gz'), 'wt', compression='gzip') as out:
        out.write(ujson.dumps(schema, double_precision=2, orient='values'))

    pg_marker_keys = get_pegasus_marker_keys(adata)
    for key in list(adata.varm.keys()):
        if key not in pg_marker_keys:
            del adata.varm[key]

    sc_marker_keys = get_scanpy_marker_keys(adata)

    for key in list(adata.uns.keys()):
        if key == 'modules':
            continue
        keep = False
        if key in sc_marker_keys:
            keep = True
        elif key.endswith('_colors'):
            field = key[0:len(key) - len('_colors')]
            if field in dataset.obs:
                keep = True
        if not keep:
            del adata.uns[key]
    adata.write(os.path.join(output_directory, 'data.h5ad'))
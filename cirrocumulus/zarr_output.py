import pandas._libs.json as ujson
import scipy.sparse
import zarr
from anndata._io.zarr import write_attribute

from cirrocumulus.anndata_util import DataType, get_scanpy_marker_keys, get_pegasus_marker_keys


def save_datasets_zarr(datasets, schema, output_directory, filesystem, whitelist):
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

    uns_whitelist = set(['modules', 'cirro-schema'])
    adata.uns['cirro-schema'] = ujson.dumps(schema, double_precision=2, orient='values')
    dataset_kwargs = {}
    chunks = None
    group = zarr.open_group(filesystem.get_mapper(output_directory), mode='a')
    if whitelist is None or 'X' in whitelist:
        if chunks is not None and not scipy.sparse.issparse(adata.X):
            write_attribute(group, "X", adata.X, dict(chunks=chunks, **dataset_kwargs))
        else:
            write_attribute(group, "X", adata.X, dataset_kwargs)
        # if scipy.sparse.issparse(adata.X):
        #     write_csc(group, 'X', adata.X)
        # else:
        #     write_attribute(group, "X", adata.X, dict(chunks=chunks, **dataset_kwargs))
        if module_dataset is not None:
            write_attribute(group, "uns/modules/X", module_dataset.X, dict(chunks=chunks, **dataset_kwargs))
            write_attribute(group, "uns/modules/var", module_dataset.var, dataset_kwargs)
    if whitelist is None or 'obs' in whitelist:
        write_attribute(group, "obs", adata.obs, dataset_kwargs)
    if whitelist is None or 'obsm' in whitelist:
        write_attribute(group, "obsm", adata.obsm, dataset_kwargs)

    pg_marker_keys = get_pegasus_marker_keys(adata)
    for key in list(adata.varm.keys()):
        if key not in pg_marker_keys:
            del adata.varm[key]
    write_attribute(group, 'varm', adata.varm, dataset_kwargs)
    write_attribute(group, "var", adata.var, dataset_kwargs)
    # keep DE results and colors
    sc_marker_keys = get_scanpy_marker_keys(adata)
    for key in list(adata.uns.keys()):
        if key in uns_whitelist:
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

    for key in list(adata.uns.keys()):
        # need to write individual groups so don't overwrite uns
        write_attribute(group, "uns/{}".format(key), adata.uns[key], dataset_kwargs)

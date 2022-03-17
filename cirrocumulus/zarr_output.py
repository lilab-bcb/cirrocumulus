import pandas._libs.json as ujson
import zarr

from cirrocumulus.anndata_util import get_pegasus_marker_keys, ADATA_MODULE_UNS_KEY
from cirrocumulus.anndata_zarr import write_attribute


def save_dataset_zarr(dataset, schema, output_directory, filesystem, whitelist):
    module_dataset = None
    if dataset.uns.get(ADATA_MODULE_UNS_KEY) is not None:
        module_dataset = dataset.uns[ADATA_MODULE_UNS_KEY]
        module_dataset.var.index.name = 'id'

    dataset.obs.index.name = 'id'
    dataset.var.index.name = 'id'
    dataset.strings_to_categoricals()
    if module_dataset is not None:
        module_dataset.strings_to_categoricals()

    dataset.uns['cirro-schema'] = ujson.dumps(schema, double_precision=2, orient='values')
    group = zarr.open_group(filesystem.get_mapper(output_directory), mode='a')

    if whitelist is None or 'X' in whitelist:
        write_attribute(group, "X", dataset.X)
        # if scipy.sparse.issparse(dataset.X):
        #     write_csc(group, 'X', dataset.X)
        # else:
        #     write_attribute(group, "X", dataset.X, dict(chunks=chunks, **dataset_kwargs))
        if module_dataset is not None:
            write_attribute(group, "uns/module/X", module_dataset.X)
            write_attribute(group, "uns/module/var", module_dataset.var)
    if whitelist is None or 'obs' in whitelist:
        write_attribute(group, "obs", dataset.obs)
    if whitelist is None or 'obsm' in whitelist:
        write_attribute(group, "obsm", dataset.obsm)

    pg_marker_keys = get_pegasus_marker_keys(dataset)
    for key in list(dataset.varm.keys()):
        if key not in pg_marker_keys:
            del dataset.varm[key]
    write_attribute(group, 'varm', dataset.varm)
    write_attribute(group, "var", dataset.var)
    # uns_whitelist = set(['module', 'cirro-schema'])
    # keep DE results and colors
    # sc_marker_keys = get_scanpy_marker_keys(dataset)
    # for key in list(dataset.uns.keys()):
    #     if key in uns_whitelist:
    #         continue
    #     keep = False
    #     if key in sc_marker_keys:
    #         keep = True
    #     elif key.endswith('_colors'):
    #         field = key[0:len(key) - len('_colors')]
    #         if field in dataset.obs:
    #             keep = True
    #     if not keep:
    #         del dataset.uns[key]

    for key in list(dataset.uns.keys()):
        # need to write individual groups so don't overwrite uns
        write_attribute(group, "uns/{}".format(key), dataset.uns[key])

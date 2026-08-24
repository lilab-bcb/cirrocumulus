import zarr
from anndata.io import write_elem

from cirrocumulus.anndata_util import ADATA_MODULE_UNS_KEY, get_pegasus_marker_keys
from cirrocumulus.util import dumps


def save_dataset_zarr(dataset, schema, output_directory, filesystem, whitelist):
    module_dataset = None
    if dataset.uns.get(ADATA_MODULE_UNS_KEY) is not None:
        module_dataset = dataset.uns[ADATA_MODULE_UNS_KEY]
        module_dataset.var.index.name = "id"

    dataset.obs.index.name = "id"
    dataset.var.index.name = "id"
    dataset.strings_to_categoricals()
    if module_dataset is not None:
        module_dataset.strings_to_categoricals()

    dataset.uns["cirro-schema"] = dumps(schema, double_precision=2, orient="values")
    # zarr resolves the protocol in the path the same way get_fs() does, and unlike an fsspec
    # mapper it creates the directories it writes into
    group = zarr.open_group(output_directory, mode="a")

    if whitelist["x"]:
        write_elem(group, "X", dataset.X)
        layers = group.require_group("layers")
        for layer in dataset.layers.keys():
            write_elem(layers, layer, dataset.layers[layer])
        if module_dataset is not None:
            module = group.require_group("uns").require_group("module")
            write_elem(module, "X", module_dataset.X)
            write_elem(module, "var", module_dataset.var)
    if whitelist["obs"]:
        write_elem(group, "obs", dataset.obs)
    if whitelist["obsm"]:
        write_elem(group, "obsm", dict(dataset.obsm))

    pg_marker_keys = get_pegasus_marker_keys(dataset)
    for key in list(dataset.varm.keys()):
        if key not in pg_marker_keys:
            del dataset.varm[key]
    write_elem(group, "varm", dict(dataset.varm))
    write_elem(group, "var", dataset.var)

    # write uns entries individually so that a rerun does not replace the whole group
    uns = group.require_group("uns")
    for key in list(dataset.uns.keys()):
        write_elem(uns, key, dataset.uns[key])

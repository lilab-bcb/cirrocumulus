import json
import os

import pandas as pd
import scipy.sparse
import tiledb
from anndata import AnnData

from cirrocumulus.abstract_dataset import AbstractDataset


class TileDBDataset(AbstractDataset):

    def get_suffixes(self):
        return ['cxg']

    def get_schema(self, filesystem, path):
        # path is to directory
        schema_dict = {'version': '1.0.0'}
        schema_dict['markers'] = []

        with tiledb.Array(os.path.join(path, 'X'), mode="r") as array:
            schema_dict['shape'] = array.shape
        annotations = {}
        for ax in ["obs"]:
            with tiledb.open(os.path.join(path, ax), mode="r") as array:
                schema_hints = json.loads(array.meta["cxg_schema"]) if "cxg_schema" in array.meta else {}
                if type(schema_hints) is not dict:
                    raise TypeError("Array schema was malformed.")
                cols = []
                for attr in array.schema:
                    schema = dict(name=attr.name, writable=False)
                    type_hint = schema_hints.get(attr.name, {})
                    # type hints take precedence
                    if "type" in type_hint:
                        schema["type"] = type_hint["type"]
                        if schema["type"] == "categorical" and "categories" in type_hint:
                            schema["categories"] = type_hint["categories"]
                    # else:
                    #     schema.update(get_schema_type_hint_from_dtype(attr.dtype))
                    cols.append(schema)

                annotations[ax] = dict(columns=cols)

                if "index" in schema_hints:
                    annotations[ax].update({"index": schema_hints["index"]})
        obs = []
        obs_cat = []
        category_order = {}
        for c in annotations['obs']['columns']:
            if c['name'] == 'name_0':  # index
                continue
            if 'type' in c and c['type'] == 'categorical':
                obs_cat.append(c['name'])
                category_order[c['name']] = c['categories']
            else:
                obs.append(c['name'])
        schema_dict['obsIndex'] = annotations['obs'].get('index', 'name_0')
        schema_dict['obs'] = obs
        schema_dict['obsCat'] = obs_cat
        schema_dict['categoryOrder'] = category_order

        # with tiledb.Array(os.path.join(path, 'cxg_group_metadata'), mode="r") as gmd:
        #     # cxg_version = gmd.meta["cxg_version"]
        #     # # version 0.1 used a malformed/shorthand semver string.
        #     # if cxg_version == "0.1" or cxg_version == "0.2.0":
        #     #     cxg_properties = json.loads(gmd.meta["cxg_properties"])
        #     colors = json.loads(gmd.meta["cxg_category_colors"]) if "cxg_category_colors" in gmd.meta else dict()
        embeddings_path_type = []
        tiledb.ls(os.path.join(path, 'emb'), lambda path, type: embeddings_path_type.append((path, type)))
        embeddings = []
        schema_dict['embeddings'] = embeddings
        for path_type in embeddings_path_type:
            if path_type[1] == 'array':
                with tiledb.open(path_type[0], mode="r") as array:
                    name = os.path.basename(path_type[0])
                    dimensions = array.ndim
                    if dimensions > 2:
                        embeddings.append({"name": name, "dimensions": 3})
                        embeddings.append({"name": name, "dimensions": 2})
                    elif dimensions == 2:
                        embeddings.append({"name": name, "dimensions": 2})

        with tiledb.open(os.path.join(path, 'var'), mode="r") as array:
            schema_dict['var'] = pd.Index(array.query(attrs=['name_0'])[:]['name_0'])
        return schema_dict

    def read_dataset(self, filesystem, path, keys=None, dataset=None):
        keys = keys.copy()
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis_keys = keys.pop('basis', [])
        dataset_info = self.get_dataset_info(filesystem, path)
        X = None
        obs = None
        var = None
        obsm = {}
        if len(var_keys) > 0:
            with tiledb.open(os.path.join(path, 'X'), mode="r") as array:
                if array.schema.sparse:
                    raise ValueError('Sparse data not supported')
                var_names = dataset_info['var']
                indices = var_names.get_indexer_for(var_keys).tolist()
                X = array.multi_index[:, indices]['']
                X = scipy.sparse.csc_matrix(X)
                var = pd.DataFrame(index=var_keys)
        if len(obs_keys) > 0:
            obs = pd.DataFrame()
            _obs_keys = []

            for key in obs_keys:
                if key == 'index':
                    _obs_keys.append(dataset_info.get('obsIndex', 'name_0'))
                else:
                    _obs_keys.append(key)
            with tiledb.open(os.path.join(path, 'obs'), mode="r") as array:
                ordered_dict = array.query(attrs=_obs_keys)[:]
            for key in ordered_dict:
                obs[key] = ordered_dict[key]

        if len(basis_keys) > 0:
            for key in basis_keys:
                with tiledb.open(os.path.join(path, 'emb', key), mode="r") as array:
                    obsm[key] = array[:]
                    if X is None:
                        X = scipy.sparse.coo_matrix(([], ([], [])), shape=(array.shape[0], 0))
        return AnnData(X=X, obs=obs, var=var, obsm=obsm)

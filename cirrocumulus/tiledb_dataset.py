import json
import os

import pandas as pd
import scipy.sparse
import tiledb


class TileDBDataset:

    def get_suffixes(self):
        return ['cxg']

    def schema(self, file_system, path):
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

    def read_dataset(self, file_system, path, keys=None, dataset=None, schema=None):
        keys = keys.copy()
        var_keys = keys.pop('X', [])
        obs_keys = keys.pop('obs', [])
        basis = keys.pop('basis', [])
        df = None
        if len(var_keys) > 0:
            with tiledb.open(os.path.join(path, 'X'), mode="r") as array:
                if array.schema.sparse:
                    raise ValueError('Sparse data not supported')
                var_names = schema['var']
                indices = var_names.get_indexer_for(var_keys).tolist()
                X = array.multi_index[:, indices]['']
                if scipy.sparse.issparse(X):
                    df = pd.DataFrame.sparse.from_spmatrix(X, columns=var_keys)
                else:  # force sparse
                    df = pd.DataFrame()
                    for i in range(len(var_keys)):
                        df[var_keys[i]] = pd.arrays.SparseArray(X[:, i], fill_value=0)

            # for key in keys.keys(): uns not supported
        #     if df is None:
        #         df = pd.DataFrame()
        #     d = adata.uns[key]
        #     features = keys[key]
        #     X = d['X'][:, d['var'].index.get_indexer_for(features)]
        #     for i in range(len(features)):
        #         df[features[i]] = X[:, i]
        if len(obs_keys) > 0:
            if df is None:
                df = pd.DataFrame()
            _obs_keys = []

            for key in obs_keys:
                if key == 'index':
                    _obs_keys.append(schema.get('obsIndex', 'name_0'))
                else:
                    _obs_keys.append(key)
            with tiledb.open(os.path.join(path, 'obs'), mode="r") as array:
                ordered_dict = array.query(attrs=_obs_keys)[:]
            for key in ordered_dict:
                df[key] = ordered_dict[key]

        if basis is not None and len(basis) > 0:
            if df is None:
                df = pd.DataFrame()
            for b in basis:
                embedding_name = b['name']
                with tiledb.open(os.path.join(path, 'emb', embedding_name), mode="r") as array:
                    dimensions = b['dimensions']
                    for i in range(dimensions):
                        df[b['coordinate_columns'][i]] = array[:, i]
        if df is None:
            df = pd.DataFrame(index=pd.RangeIndex(schema['shape'][0]))
        return df

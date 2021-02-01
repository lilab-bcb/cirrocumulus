import argparse
import gzip
import json
import logging
import os

import numpy as np
import pandas as pd
import pandas._libs.json as ujson
import scipy.sparse
from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.dataset_api import DatasetAPI
from cirrocumulus.io_util import get_markers, filter_markers, add_spatial, SPATIAL_HELP, unique_id
from cirrocumulus.simple_data import SimpleData

logger = logging.getLogger("cirro")

cluster_fields = ['seurat_clusters', 'leiden', 'louvain']
categorical_fields_convert = ['seurat_clusters']


def read_adata(path, backed=False, spatial_directory=None, use_raw=False):
    import anndata
    adata = anndata.read_loom(path) if path.lower().endswith('.loom') else anndata.read(path,
        backed=backed)
    if use_raw and adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
        logger.info('Using adata.raw')
        adata = anndata.AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns)

    if spatial_directory is not None:
        if not add_spatial(adata, spatial_directory):
            print('No spatial data found in {}'.format(spatial_directory))
    if not backed:
        if scipy.sparse.issparse(adata.X) and scipy.sparse.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsc()

    def fix_column_names(df):
        rename = {}
        for c in df.columns:
            if c.find(' ') != -1:
                rename[c] = c.replace(' ', '_')
        return df.rename(rename, axis=1) if len(rename) > 0 else df

    adata.obs = fix_column_names(adata.obs)
    adata.var = fix_column_names(adata.var)
    for field in categorical_fields_convert:
        if field in adata.obs and not pd.api.types.is_categorical_dtype(adata.obs[field]):
            logger.info('Converting {} to categorical'.format(field))
            adata.obs[field] = adata.obs[field].astype('category')
    for key in adata.obsm:
        if key.find(' ') != -1:
            new_key = key.replace(' ', '_')
            adata.obsm[new_key] = adata.obsm[key]
            del adata.obsm[key]
    return adata


def to_json(data):
    return ujson.dumps(data, double_precision=2, orient='values')


def make_unique(index, join='-1'):
    index = index.str.replace('/', '_')
    lower_index = index.str.lower()
    if lower_index.is_unique:
        return index
    from collections import defaultdict

    indices_dup = lower_index.duplicated(keep="first")
    values_dup_lc = lower_index.values[indices_dup]
    values_dup = index.values[indices_dup]
    counter = defaultdict(lambda: 0)
    for i, v in enumerate(values_dup_lc):
        counter[v] += 1
        values_dup[i] += join + str(counter[v])
    values = index.values
    values[indices_dup] = values_dup
    index = pd.Index(values)
    return index


def write_obs_stats(directory, results):
    os.makedirs(os.path.join(directory, 'stats', 'obs'), exist_ok=True)
    for column in results:
        dimension_or_measure_summary = results[column]
        if 'categories' in dimension_or_measure_summary:
            with gzip.open(os.path.join(directory, 'stats', 'obs', column + '.json.gz'), 'wt') as f:
                f.write(to_json(dimension_or_measure_summary))
        else:
            with gzip.open(os.path.join(directory, 'stats', 'obs', column + '.json.gz'), 'wt') as f:
                f.write(to_json(dimension_or_measure_summary))


def write_X_bins(base_output_dir, n_bins=25):
    # read in means from stats/X
    obs_avg = []
    names = []
    d = os.path.join(base_output_dir, 'stats', 'X')
    for file in os.listdir(d):
        if file.endswith('.json.gz'):
            with gzip.open(os.path.join(d, file), 'rt') as f:
                names.append(file[0:file.rindex('.json.gz')])
                result = json.load(f)
                obs_avg.append(result['mean'])

    obs_avg = pd.Series(np.array(obs_avg), index=names)
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    obs_cut = obs_cut.astype(int)
    with gzip.open(os.path.join(base_output_dir, 'X.json.gz'), 'wt') as f:
        f.write(to_json(obs_cut.to_dict()))


def write_X_stats(directory, results):
    os.makedirs(os.path.join(directory, 'stats', 'X'), exist_ok=True)
    for column in results:
        measure_summary = results[column]
        with gzip.open(os.path.join(directory, 'stats', 'X', column + '.json.gz'), 'wt') as f:
            f.write(to_json(measure_summary))


def write_grouped_stats(directory, dotplot_results):
    # mean and fraction expressed for every gene per category
    for dotplot_result in dotplot_results:
        categories = dotplot_result['categories']
        category_name = dotplot_result['name']
        os.makedirs(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'X'), exist_ok=True)
        with gzip.open(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'index.json.gz'),
                'wt') as f:
            f.write(to_json(categories))
        values = dotplot_result['values']
        for value in values:
            var_field = value['name']
            with gzip.open(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'X', var_field + '.json.gz'),
                    'wt') as f:
                f.write(to_json(value))


def write_basis_X(coords_group, basis_group, adata, result):
    shape = (coords_group['index'].shape[0], adata.shape[1])
    X = basis_group.require_dataset('X', shape=shape, chunks=(None, 1), dtype='f4')
    for column in result['values']:
        X_index = adata.var.index.get_indexer_for([column])[0]
        X[:, X_index] = result['values'][column]


def require_binned_basis_group(store, basis):
    basis_group = store.require_group('obsm_summary/' + basis['full_name'])
    basis_group.attrs['nbins'] = basis['nbins']
    basis_group.attrs['agg'] = basis['agg']
    basis_group.attrs['name'] = basis['name']
    basis_group.attrs['dimensions'] = basis['dimensions']
    basis_group.attrs['precomputed'] = True
    return basis_group


def write_basis_obs(basis, coords_group, obs_group, result):
    coords_group['index'] = result['bins']
    ndim = len(basis['coordinate_columns'])
    nbins = len(result['bins'])

    value_ds = coords_group.require_dataset('value', dtype='f4', shape=(nbins, ndim))

    for i in range(len(basis['coordinate_columns'])):
        coord_field = basis['coordinate_columns'][i]
        value_ds[:, i] = np.array(result['coordinates'][coord_field])
    # write_pq(bin_dict, output_directory, 'index')

    for column in result['values']:
        vals = result['values'][column]
        if not isinstance(vals, dict):  # continuous field
            obs_group[column] = vals
        else:
            g = obs_group.require_group(column)
            g['purity'] = vals['purity']
            g['mode'] = vals['value']
        # write_pq(vals, output_directory, column)


class PrepareData:

    def __init__(self, adata, output, dimensions=None, groups=[], group_nfeatures=10, markers=[],
                 output_format='parquet'):
        self.adata = adata
        self.groups = groups
        self.group_nfeatures = group_nfeatures
        self.markers = markers
        self.stats = False
        self.nbins = None
        self.bin_agg_function = None
        self.output_format = output_format

        # if basis_list is None or len(basis_list) == 0:
        #     basis_list = list(self.adata.obsm_keys())
        self.basis_list_to_precompute = None

        index = make_unique(self.adata.var.index.append(pd.Index(self.adata.obs.columns)))
        self.adata.var.index = index[0:len(self.adata.var.index)]
        self.adata.obs.columns = index[len(self.adata.var.index):]

        self.base_output_dir = output

        self.dataset_api = DatasetAPI()
        ds = AnndataDataset()
        ds.path_to_data[''] = self.adata
        self.dataset_api.add(ds)

        self.input_dataset = {'name': '', 'url': '', 'id': ''}
        dimensions_supplied = dimensions is not None and len(dimensions) > 0
        self.dimensions = [] if not dimensions_supplied else dimensions
        self.measures = []
        self.others = []
        for i in range(len(self.adata.obs.columns)):
            name = self.adata.obs.columns[i]
            c = self.adata.obs[name]
            if pd.api.types.is_object_dtype(c):
                self.adata.obs[name] = self.adata.obs[name].astype('category')
                c = self.adata.obs[name]
            if not dimensions_supplied and pd.api.types.is_categorical_dtype(c):
                if 1 < len(c.cat.categories) < 2000:
                    self.dimensions.append(name)
                    if c.isna().sum() > 0:
                        logger.info('Replacing nans in {}'.format(name))
                        self.adata.obs[name] = self.adata.obs[name].astype(str)
                        self.adata.obs.loc[self.adata.obs[name].isna(), name] = ''
                        self.adata.obs[name] = self.adata.obs[name].astype('category')
                else:
                    self.others.append(name)
            elif not pd.api.types.is_string_dtype(c) and not pd.api.types.is_object_dtype(c):
                self.measures.append('obs/' + name)
            else:
                self.others.append(name)

    def get_path(self, path):
        return os.path.join(self.base_output_dir, path)

    def execute(self):
        basis_list = self.basis_list_to_precompute
        nbins = self.nbins
        group_nfeatures = self.group_nfeatures
        bin_agg_function = self.bin_agg_function
        output_format = self.output_format
        if not os.path.exists(self.base_output_dir):
            os.makedirs(self.base_output_dir, exist_ok=True)

        schema = self.get_schema()
        results = schema.get('results', [])
        if len(results) > 0:
            uns_dir = os.path.join(self.base_output_dir, 'uns')
            os.makedirs(uns_dir, exist_ok=True)
            for i in range(len(results)):  # keep id, name, type in schema, store rest in file
                result = results[i]
                result_id = result.pop('id')
                results[i] = dict(id=result_id, name=result.pop('name'), type=result.pop('type'),
                    content_type='application/json', content_encoding='gzip')
                with gzip.open(os.path.join(uns_dir, result_id + '.json.gz'), 'wt') as f:
                    f.write(to_json(result))

        markers = schema.get('markers', [])

        if len(markers) == 0 and self.groups is None:
            groups = []
            for field in self.adata.obs.columns:
                field_lc = field.lower()
                for cluster_field in cluster_fields:
                    if field_lc.find(cluster_field) != -1:
                        groups.append(field)
                        break
            self.groups = groups
        if self.groups is not None:
            for field in self.groups:
                if not pd.api.types.is_categorical_dtype(self.adata.obs[field]):
                    self.adata.obs[field] = self.adata.obs[field].astype('category')
                if len(self.adata.obs[field].cat.categories) > 1:
                    markers += SimpleData.find_markers(self.adata, field, group_nfeatures)
            schema['markers'] = markers

        images = self.adata.uns.get('images')
        if images is not None:
            image_dir = os.path.join(self.base_output_dir, 'images')
            if not os.path.exists(image_dir):
                os.mkdir(image_dir)
            for image in images:
                path = image['image']
                import shutil
                shutil.copy(path, os.path.join(image_dir, os.path.basename(path)))
                image['image'] = 'images/' + os.path.basename(path)

        if output_format == 'parquet':
            from cirrocumulus.parquet_io import save_adata_pq
            save_adata_pq(self.adata, schema, self.base_output_dir)
        elif output_format == 'json':
            from cirrocumulus.json_io import save_adata_json
            save_adata_json(self.adata, schema, self.base_output_dir)
        elif output_format == 'jsonl':
            from cirrocumulus.jsonl_io import save_adata_jsonl
            save_adata_jsonl(self.adata, schema, self.base_output_dir)
        else:
            raise ValueError("Unknown format")

        # if self.stats:
        #     self.summary_stats()
        #     self.grouped_stats()

        if nbins is not None and basis_list is not None:
            for basis_name in basis_list:
                self.grid_embedding(basis_name, bin_agg_function, nbins)


    # def summary_stats(self):
    #     dimensions = self.dimensions
    #     measures = self.measures
    #     dataset_api = self.dataset_api
    #     input_dataset = self.input_dataset
    #     X_range = self.X_range
    #     adata = self.adata
    #     base_dir = self.base_output_dir
    #
    #     if X_range[0] == 0:
    #         process_results = data_processing.handle_data(dataset_api=dataset_api, dataset=input_dataset,
    #             stats=dict(measures=measures, dimensions=dimensions))
    #         write_obs_stats(base_dir, process_results['summary'])
    #
    #     logger.info('Summary stats X {}-{}'.format(X_range[0], X_range[1]))
    #     process_results = data_processing.handle_data(dataset_api=dataset_api, dataset=input_dataset,
    #         stats=dict(measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[]))
    #     write_X_stats(base_dir, process_results['summary'])
    #     write_X_bins(base_dir)
    #
    # # stats, grouped by categories in adata.obs for dot plot
    # def grouped_stats(self):
    #     dimensions = self.dimensions
    #     dataset_api = self.dataset_api
    #     input_dataset = self.input_dataset
    #     adata = self.adata
    #     X_range = self.X_range
    #     base_dir = self.base_output_dir
    #     logger.info('Grouped stats X {}-{}'.format(X_range[0], X_range[1]))
    #     process_results = data_processing.handle_data(dataset_api=dataset_api, dataset=input_dataset,
    #         grouped_stats=dict(measures=adata.var_names[X_range[0]:X_range[1]], dimensions=dimensions))
    #     dotplot_results = process_results['dotplot']
    #
    #     write_grouped_stats(base_dir, dotplot_results)


    # def grid_embedding(self, basis_name, summary, nbins):
    #     if nbins <= 0:
    #         return
    #     logger.info('{} embedding'.format(basis_name))
    #     basis = get_basis(basis_name, nbins, summary, 2, False)
    #
    #     dimensions = self.dimensions
    #     measures = self.measures
    #     dataset_api = self.dataset_api
    #     input_dataset = self.input_dataset
    #     adata = self.adata
    #     X_range = self.X_range
    #     full_basis_name = basis['full_name']
    #
    #     basis_group = require_binned_basis_group(store, basis)
    #     obs_group = basis_group.require_group('obs')
    #     coords_group = basis_group.require_group('coords')
    #
    #     if X_range[0] == 0:
    #         df_with_coords = pd.DataFrame()
    #         obsm = adata.obsm[basis_name]
    #         for i in range(len(basis['coordinate_columns'])):
    #             df_with_coords[basis['coordinate_columns'][i]] = obsm[:, i]
    #         EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'],
    #             bin_name=full_basis_name)
    #         coords_group['bins'] = df_with_coords[full_basis_name].values
    #         coords_group['bin_coords'] = df_with_coords[basis['coordinate_columns']].values
    #
    #         result = data_processing.handle_data(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
    #             measures=measures + ['__count'], dimensions=dimensions, quick=False)
    #         write_basis_obs(basis, coords_group, obs_group, result)
    #         logger.info('{} embedding finished writing obs'.format(basis_name))
    #     # write X to obsm_summary/name
    #     logger.info('{} embedding X {}-{}'.format(basis_name, X_range[0], X_range[1]))
    #     result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
    #         measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[])
    #     write_basis_X(coords_group, basis_group, adata, result)
    #     # write_pq(dict(value=result['values'][column]), output_directory, column)

    def get_schema(self):
        basis_list = self.basis_list_to_precompute
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        result = SimpleData.schema(self.adata)
        markers = result.get('markers', [])

        if self.markers is not None:  # add from file
            markers += get_markers(self.markers)

        markers = filter_markers(self.adata, markers)

        for marker in markers:
            if marker.get('id') is None:
                marker['id'] = unique_id()
            marker['readonly'] = True
        result['markers'] = markers
        result['format'] = self.output_format
        if self.stats:
            result['precomputed'] = True  # has stats
        if nbins is not None:
            result['embeddings'] = []
            for basis_name in basis_list:
                result['embeddings'].append(
                    {'precomputed': True, 'name': basis_name, 'nbins': nbins, 'agg': bin_agg_function,
                     'dimensions': 2})
        return result


def main(argsv):
    parser = argparse.ArgumentParser(
        description='Prepare a dataset for cirrocumulus server.')
    parser.add_argument('dataset', help='Path to a h5ad, loom, or Seurat file')
    parser.add_argument('--out', help='Path to output directory')
    # parser.add_argument('--stats', dest="stats", help='Generate precomputed stats', action='store_true')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    # parser.add_argument('--basis', help='List of embeddings to precompute', action='append')
    #  parser.add_argument('--groups', help='List of groups to precompute summary statistics', action='append')
    parser.add_argument('--markers',
        help='Path to JSON file of precomputed markers that maps name to features. For example {"a":["gene1", "gene2"], "b":["gene3"]',
        action='append')
    parser.add_argument('--groups',
        help='List of groups to compute markers for (e.g. louvain). Note that markers created with scanpy or cumulus are automatically included.',
        action='append')
    parser.add_argument('--group_nfeatures', help='Number of marker genes/features to include', type=int, default=10)
    parser.add_argument('--spatial', help=SPATIAL_HELP)
    # parser.add_argument('--output_format', help='Output file format', choices=['json', 'parquet'], default='parquet')
    # parser.add_argument('--nbins', help='Number of bins. Set to 0 to disable binning', default=500, type=int)
    # parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    # parser.add_argument('--X_range', help='Start and end position of data matrix (e.g. 0-5000)', type=str)

    args = parser.parse_args(argsv)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    input_basis = []  # check if specified as comma delimited list
    input_X_range = None
    out = args.out
    input_dataset = args.dataset
    output_format = 'parquet'  # args.output_format
    if out is None:
        out = os.path.basename(input_dataset)
        out = out[0:out.rindex('.')]
    tmp_file = None
    use_raw = False
    if input_dataset.lower().endswith('.rds'):
        import subprocess
        import tempfile
        import pkg_resources

        _, h5_file = tempfile.mkstemp(suffix='.h5ad')
        os.remove(h5_file)
        subprocess.check_call(
            ['Rscript', pkg_resources.resource_filename("cirrocumulus", 'seurat2h5ad.R'), input_dataset, h5_file])
        input_dataset = h5_file
        tmp_file = h5_file
        use_raw = True
    # summary = 'max'
    # nbins = -1
    # summary = args.summary
    # if args.basis is not None:
    #     for b in args.basis:
    #         input_basis.extend(b.split(','))
    # input_X_range = args.X_range
    # if input_X_range is not None:
    #     input_X_range = input_X_range.split('-')
    #     assert len(input_X_range) == 2
    #     input_X_range[0] = int(input_X_range[0])
    #     input_X_range[1] = int(input_X_range[1])

    adata = read_adata(input_dataset, backed=args.backed, spatial_directory=args.spatial, use_raw=use_raw)
    prepare_data = PrepareData(adata=adata, output=out, dimensions=args.groups, groups=args.groups,
        group_nfeatures=args.group_nfeatures,
        markers=args.markers, output_format=output_format)
    prepare_data.execute()
    if tmp_file is not None:
        os.remove(tmp_file)


if __name__ == '__main__':
    import sys

    main(sys.argv)

import argparse
import gzip
import json
import logging
import os

import anndata
import numpy as np
import pandas as pd
import zarr
from natsort import natsorted
from pandas import CategoricalDtype

import cirrocumulus.data_processing as data_processing
from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.dataset_api import DatasetAPI
from cirrocumulus.embedding_aggregator import EmbeddingAggregator, get_basis
from cirrocumulus.entity import Entity
from cirrocumulus.simple_data import SimpleData

logger = logging.getLogger("cirro")


def make_unique(index, join='-1'):
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


def make_ordered(df, columns):
    if columns is None:
        columns = df.columns
    for i in range(len(columns)):
        name = columns[i]
        c = df[name]
        if pd.api.types.is_categorical_dtype(c) and not c.dtype.ordered:
            df[name] = df[name].astype(CategoricalDtype(natsorted(df[name].dtype.categories), ordered=True))


def write_obs_stats(directory, results):
    os.makedirs(os.path.join(directory, 'stats', 'obs'), exist_ok=True)
    for column in results:
        dimension_or_measure_summary = results[column]
        if 'categories' in dimension_or_measure_summary:
            with gzip.open(os.path.join(directory, 'stats', 'obs', column + '.json.gz'), 'wt') as f:
                json.dump(dimension_or_measure_summary, f)
            # write_pq(dict(index=dimension_or_measure_summary['categories'],
            #     value=dimension_or_measure_summary['counts']),
            #     dimension_output_directory, column)
        else:
            with gzip.open(os.path.join(directory, 'stats', 'obs', column + '.json.gz'), 'wt') as f:
                json.dump(dimension_or_measure_summary, f)


def write_X_stats(directory, results):
    os.makedirs(os.path.join(directory, 'stats', 'X'), exist_ok=True)
    for column in results:
        measure_summary = results[column]
        with gzip.open(os.path.join(directory, 'stats', 'X', column + '.json.gz'), 'wt') as f:
            json.dump(measure_summary, f)


def write_grouped_stats(directory, adata, dotplot_results):
    # mean and fraction expressed for every gene per category
    for dotplot_result in dotplot_results:
        # output_directory = self.get_path('grouped_statistics/' + dotplot_result['name'])
        # write_pq(dict(index=dotplot_result['categories']), output_directory, 'index')
        categories = dotplot_result['categories']
        category_name = dotplot_result['name']
        os.makedirs(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'X'), exist_ok=True)
        with gzip.open(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'index.json.gz'),
                'wt') as f:
            json.dump(categories, f)
        values = dotplot_result['values']
        for value in values:
            var_field = value['name']
            with gzip.open(os.path.join(directory, 'grouped_stats', 'obs', category_name, 'X', var_field + '.json.gz'),
                    'wt') as f:
                json.dump(value, f)


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

    def __init__(self, input_path, backed, basis_list, nbins, bin_agg_function, output,
                 X_range=None, dimensions=None):
        self.input_path = input_path
        self.adata = anndata.read(input_path, backed=backed)
        if not backed:
            self.adata = self.adata[:, self.adata.X.sum(axis=0).A1 > 0]
        index = make_unique(self.adata.var.index.append(pd.Index(self.adata.obs.columns)))
        self.adata.var.index = index[0:len(self.adata.var.index)]
        self.adata.obs.columns = index[len(self.adata.var.index):]

        logger.info('{} - {} x {}'.format(input_path, self.adata.shape[0], self.adata.shape[1]))
        self.X_range = (0, self.adata.shape[1]) if X_range is None else X_range
        self.base_name = output
        if basis_list is None or len(basis_list) == 0:
            basis_list = list(self.adata.obsm_keys())
        self.basis_list = basis_list
        self.nbins = nbins
        self.dataset_api = DatasetAPI()
        self.bin_agg_function = bin_agg_function
        ds = AnndataDataset()
        ds.path_to_data[input_path] = self.adata
        self.dataset_api.add(ds)

        self.input_dataset = Entity(input_path,
            {'name': os.path.splitext(os.path.basename(input_path))[0], 'url': input_path})
        dimensions_supplied = dimensions is not None and len(dimensions) > 0
        self.dimensions = [] if not dimensions_supplied else dimensions
        self.measures = []
        self.others = []
        for i in range(len(self.adata.obs.columns)):
            name = self.adata.obs.columns[i]
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
        make_ordered(self.adata.obs, self.dimensions)
        self.store = zarr.open(self.base_name)

    def get_path(self, path):
        return os.path.join(self.base_name, path)

    def execute(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        X_range = self.X_range
        if X_range[0] == 0:
            import scipy.sparse
            if scipy.sparse.issparse(self.adata.X) and scipy.sparse.isspmatrix_csr(self.adata.X):
                self.adata.X = self.adata.X.tocsc()
            self.adata.write_zarr(self.base_name, chunks=(None, 1))
        self.summary_stats()
        self.grouped_stats()

        for basis_name in basis_list:
            self.grid_embedding(basis_name, bin_agg_function, nbins)
        # self.schema()

    def summary_stats(self):
        dimensions = self.dimensions
        measures = self.measures
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        X_range = self.X_range
        adata = self.adata
        base_dir = self.base_name

        if X_range[0] == 0:
            process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
                measures=measures, dimensions=dimensions)
            write_obs_stats(base_dir, process_results['summary'])
            # write_pq(
            #     dict(min=[dimension_or_measure_summary['min']], max=[dimension_or_measure_summary['max']],
            #         sum=[dimension_or_measure_summary['sum']], mean=[dimension_or_measure_summary['mean']]),
            #     measure_output_directory, column)

        logger.info('Summary stats X {}-{}'.format(X_range[0], X_range[1]))
        process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[])
        write_X_stats(base_dir, process_results['summary'])


        # write_pq(dict(min=[measure_summary['min']],
        #     max=[measure_summary['max']],
        #     sum=[measure_summary['sum']],
        #     numExpressed=[measure_summary['numExpressed']],
        #     mean=[measure_summary['mean']]),
        #     measure_output_directory, column)


    # stats, grouped by categories in adata.obs for dot plot
    def grouped_stats(self):
        dimensions = self.dimensions
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        X_range = self.X_range
        base_dir = self.base_name
        logger.info('Grouped stats X {}-{}'.format(X_range[0], X_range[1]))
        process_results = data_processing.handle_grouped_stats(dataset_api=dataset_api, dataset=input_dataset,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=dimensions)
        dotplot_results = process_results['dotplot']

        write_grouped_stats(base_dir, adata, dotplot_results)
        # for value in values:
        # write_pq(dict(mean=value['mean'],
        #     fractionExpressed=value['fractionExpressed']),
        #     output_directory, value['name'])

    def grid_embedding(self, basis_name, summary, nbins):
        if nbins <= 0:
            return
        logger.info('{} embedding'.format(basis_name))
        basis = get_basis(basis_name, nbins, summary, 2, False)

        dimensions = self.dimensions
        measures = self.measures
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        X_range = self.X_range
        full_basis_name = basis['full_name']

        store = self.store
        basis_group = require_binned_basis_group(store, basis)
        obs_group = basis_group.require_group('obs')
        coords_group = basis_group.require_group('coords')

        if X_range[0] == 0:
            df_with_coords = pd.DataFrame()
            obsm = adata.obsm[basis_name]
            for i in range(len(basis['coordinate_columns'])):
                df_with_coords[basis['coordinate_columns'][i]] = obsm[:, i]
            EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'],
                bin_name=full_basis_name)
            coords_group['bins'] = df_with_coords[full_basis_name].values
            coords_group['bin_coords'] = df_with_coords[basis['coordinate_columns']].values

            result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
                measures=measures + ['__count'], dimensions=dimensions, quick=False)
            write_basis_obs(basis, coords_group, obs_group, result)
            logger.info('{} embedding finished writing obs'.format(basis_name))
        # write X to obsm_summary/name
        logger.info('{} embedding X {}-{}'.format(basis_name, X_range[0], X_range[1]))
        result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[])
        write_basis_X(coords_group, basis_group, adata, result)
        # write_pq(dict(value=result['values'][column]), output_directory, column)

    def schema(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        X_range = self.X_range
        if X_range[0] == 0:
            output_file = os.path.join(self.base_name, 'index.json')
            result = SimpleData.schema(self.adata)
            result['precomputed'] = True
            if nbins > 0:
                result['embeddings'] = []
                for basis_name in basis_list:
                    result['embeddings'].append(
                        {'precomputed': True, 'name': basis_name, 'nbins': nbins, 'agg': bin_agg_function,
                         'dimensions': 2})
            with gzip.open(output_file, 'wt') as f:
                json.dump(result, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='(BETA!) Prepare a dataset by binning on a grid using an embedding, generating feature statistics feature statistics within a category, and saving the data for easy slicing by feature.')
    parser.add_argument('dataset', help='Path to a h5ad file')
    parser.add_argument('out', help='Path to output directory')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--basis', help='List of embeddings to precompute', action='append')
    parser.add_argument('--groups', help='List of groups to precompute summary statistics', action='append')
    parser.add_argument('--nbins', help='Number of bins. Set to 0 to disable binning', default=500, type=int)
    parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    parser.add_argument('--X_range', help='Start and end position of data matrix (e.g. 0-5000)', type=str)

    args = parser.parse_args()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())
    logger.info('preparing ' + args.dataset + '...')
    input_basis = []  # check if specified as comma delimited list
    if args.basis is not None:
        for b in args.basis:
            input_basis.extend(b.split(','))
    input_X_range = args.X_range
    if input_X_range is not None:
        input_X_range = input_X_range.split('-')
        assert len(input_X_range) == 2
        input_X_range[0] = int(input_X_range[0])
        input_X_range[1] = int(input_X_range[1])
    prepare_data = PrepareData(args.dataset, args.backed, input_basis, args.nbins, args.summary, args.out,
        X_range=input_X_range, dimensions=args.groups)
    prepare_data.execute()

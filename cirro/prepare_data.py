import argparse
import json
import logging
import os

import anndata
import pandas as pd

import cirro.data_processing as data_processing
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import EmbeddingAggregator, get_basis
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset
from cirro.io import save_adata
from cirro.parquet_dataset import write_pq
from cirro.simple_data import SimpleData

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
        ds = H5ADDataset()
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

    def get_path(self, path):
        return os.path.join(self.base_name, path)

    def execute(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        X_range = self.X_range
        self.summary_stats()
        self.grouped_stats()
        save_adata(self.adata, self.get_path('data'), X_range=X_range)
        for basis_name in basis_list:
            self.grid_embedding(basis_name, bin_agg_function, nbins)
        self.schema()

    def summary_stats(self):
        dimensions = self.dimensions
        measures = self.measures
        dimension_output_directory = self.get_path('counts')
        measure_output_directory = self.get_path('statistics')
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        X_range = self.X_range
        adata = self.adata

        if X_range[0] == 0:
            process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
                measures=measures, dimensions=dimensions)
            summary = process_results['summary']
            for column in summary:

                dimension_or_measure_summary = summary[column]
                if 'categories' in dimension_or_measure_summary:
                    write_pq(dict(index=dimension_or_measure_summary['categories'],
                        value=dimension_or_measure_summary['counts']),
                        dimension_output_directory, column)
                else:

                    write_pq(
                        dict(min=[dimension_or_measure_summary['min']], max=[dimension_or_measure_summary['max']],
                            sum=[dimension_or_measure_summary['sum']], mean=[dimension_or_measure_summary['mean']]),
                        measure_output_directory, column)

        logger.info('Summary stats X {}-{}'.format(X_range[0], X_range[1]))
        process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[])
        summary = process_results['summary']
        for column in summary:
            measure_summary = summary[column]
            write_pq(dict(min=[measure_summary['min']],
                max=[measure_summary['max']],
                sum=[measure_summary['sum']],
                numExpressed=[measure_summary['numExpressed']],
                mean=[measure_summary['mean']]),
                measure_output_directory, column)


    # stats, grouped by categories in adata.obs for dot plot
    def grouped_stats(self):
        dimensions = self.dimensions
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        X_range = self.X_range

        logger.info('Grouped stats X {}-{}'.format(X_range[0], X_range[1]))
        process_results = data_processing.handle_grouped_stats(dataset_api=dataset_api, dataset=input_dataset,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=dimensions)
        dotplot_results = process_results['dotplot']

        for dotplot_result in dotplot_results:
            output_directory = self.get_path('grouped_statistics/' + dotplot_result['name'])
            write_pq(dict(index=dotplot_result['categories']), output_directory, 'index')
            values = dotplot_result['values']
            for value in values:
                write_pq(dict(mean=value['mean'],
                    fractionExpressed=value['fractionExpressed']),
                    output_directory, value['name'])

    def grid_embedding(self, basis_name, summary, nbins):
        logger.info('{} embedding'.format(basis_name))
        basis = get_basis(basis_name, nbins, summary, 2, False)

        dimensions = self.dimensions
        measures = self.measures
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        X_range = self.X_range
        full_basis_name = basis['full_name']
        output_directory = self.get_path('obsm_summary/' + full_basis_name)

        # write cell level bin to /data
        if X_range[0] == 0:
            df_with_coords = pd.DataFrame()
            obsm = adata.obsm[basis_name]
            for i in range(len(basis['coordinate_columns'])):
                df_with_coords[basis['coordinate_columns'][i]] = obsm[:, i]
            EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'],
                bin_name=full_basis_name)
            dict_with_coords = df_with_coords.to_dict(orient='list')
            write_pq(dict_with_coords, os.path.join(self.base_name, 'data'), full_basis_name)
            if nbins <= 0:
                return
            # write bins and coordinates to obsm_summary/name
            result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
                measures=measures + ['__count'], dimensions=dimensions)

            bin_dict = dict(index=result['bins'])
            for column in result['coordinates']:
                bin_dict[column] = result['coordinates'][column]
            write_pq(bin_dict, output_directory, 'index')

            for column in result['values']:
                vals = result['values'][column]
                if not isinstance(vals, dict):
                    vals = dict(value=vals)
                write_pq(vals, output_directory, column)
            logger.info('{} embedding finished writing obs'.format(basis_name))
        # write X to obsm_summary/name
        logger.info('{} embedding X {}-{}'.format(basis_name, X_range[0], X_range[1]))
        result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
            measures=adata.var_names[X_range[0]:X_range[1]], dimensions=[])
        for column in result['values']:
            write_pq(dict(value=result['values'][column]), output_directory, column)

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
            with open(output_file, 'wt') as f:
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

import argparse
import json
import logging
import os

import anndata
import pandas as pd

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import EmbeddingAggregator, get_basis
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset
from cirro.io import write_table, save_adata
from cirro.simple_data import SimpleData

logger = logging.getLogger("cirro")


class PrepareData:

    def __init__(self, input_path, backed, basis_list, nbins, bin_agg_function, column_batch_size=1000):
        self.input_path = input_path
        self.column_batch_size = column_batch_size
        self.adata = anndata.read(input_path, backed=backed)
        dimensions = []
        measures = []
        others = []
        if basis_list is None:
            basis_list = list(self.adata.obsm_keys())
            for key in basis_list:
                if self.adata.obsm[key].shape[1] >= 10:
                    basis_list.remove(key)
        self.basis_list = basis_list
        self.nbins = nbins
        self.dataset_api = DatasetAPI()
        self.bin_agg_function = bin_agg_function
        ds = H5ADDataset()
        ds.path_to_data[input_path] = self.adata
        self.dataset_api.add(ds)

        self.input_dataset = Entity(input_path,
            {'name': os.path.splitext(os.path.basename(input_path))[0], 'url': input_path})
        for i in range(len(self.adata.obs.columns)):
            name = self.adata.obs.columns[i]
            c = self.adata.obs[name]
            if pd.api.types.is_categorical_dtype(c):
                if 1 < len(c.cat.categories) < 5000:
                    dimensions.append(name)
                else:
                    others.append(name)
            elif not pd.api.types.is_string_dtype(c) and not pd.api.types.is_object_dtype(c):
                measures.append('obs/' + name)
            else:
                others.append(name)
        self.others = others
        self.dimensions = dimensions
        self.measures = measures
        self.base_name = os.path.splitext(os.path.basename(input_path))[0]

    def execute(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        save_adata(self.adata, os.path.join(self.base_name, 'data'))
        self.summary_stats()
        self.grouped_stats()
        for basis_name in basis_list:
            self.grid_embedding(basis_name, bin_agg_function, nbins)
        self.schema()

    def summary_stats(self):
        dimensions = self.dimensions
        measures = self.measures
        dimension_output_directory = os.path.join(self.base_name, 'counts')
        measure_output_directory = os.path.join(self.base_name, 'statistics')
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        column_batch_size = self.column_batch_size

        process_results = process_data(dataset_api=dataset_api, dataset=input_dataset, summary_measures=measures,
            summary_dimensions=dimensions, return_types=['summary'])
        summary = process_results['summary']
        for column in summary:
            dimension_or_measure_summary = summary[column]
            if 'categories' in dimension_or_measure_summary:
                write_table(dict(index=dimension_or_measure_summary['categories'],
                    value=dimension_or_measure_summary['counts']),
                    dimension_output_directory, column, write_statistics=False)
            else:
                write_table(dict(min=[dimension_or_measure_summary['min']], max=[dimension_or_measure_summary['max']],
                    sum=[dimension_or_measure_summary['sum']], mean=[dimension_or_measure_summary['mean']]),
                    measure_output_directory, column, write_statistics=False)
        for i in range(0, self.adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            logger.info('Summary stats {}-{}'.format(i, end))
            process_results = process_data(dataset_api=dataset_api, dataset=input_dataset,
                summary_measures=adata.var_names[i:end],
                summary_dimensions=[], return_types=['summary'])
            summary = process_results['summary']
            for column in summary:
                measure_summary = summary[column]
                write_table(dict(min=[measure_summary['min']], max=[measure_summary['max']],
                    sum=[measure_summary['sum']],
                    numExpressed=[measure_summary['numExpressed']],
                    mean=[measure_summary['mean']]),
                    measure_output_directory, column, write_statistics=False)


    # stats, grouped by categories in adata.obs for dot plot
    def grouped_stats(self):
        dimensions = self.dimensions
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        column_batch_size = self.column_batch_size

        for i in range(0, self.adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            logger.info('Grouped stats {}-{}'.format(i, end))
            process_results = process_data(dataset_api=dataset_api, dataset=input_dataset,
                summary_measures=adata.var_names[i:end],
                summary_dimensions=dimensions, return_types=['dotplot'])
            dotplot_results = process_results['dotplot']

            for dotplot_result in dotplot_results:
                output_directory = os.path.join(self.base_name, 'grouped_statistics', dotplot_result['name'])
                write_table(dict(index=dotplot_result['categories']),
                    output_directory, 'index')
                values = dotplot_result['values']
                for value in values:
                    write_table(dict(mean=[value['mean']],
                        fractionExpressed=[value['fractionExpressed']]),
                        output_directory, value['name'])

    def grid_embedding(self, basis_name, summary, nbins):
        logger.info('{} embedding'.format(basis_name))
        basis = get_basis(basis_name)
        path = basis_name + '_' + str(nbins) + '_' + str(summary)
        dimensions = self.dimensions
        measures = self.measures
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        column_batch_size = self.column_batch_size
        output_directory = os.path.join(self.base_name, 'obsm_summary', path)

        # write cell level bin to /data
        df_with_coords = pd.DataFrame()
        obsm = adata.obsm[basis_name]
        for i in range(len(basis['coordinate_columns'])):
            df_with_coords[basis['coordinate_columns'][i]] = obsm[:, i]
        EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'], None)
        dict_with_coords = df_with_coords.to_dict(orient='list')
        dict_with_coords['index'] = df_with_coords.index.values
        write_table(dict_with_coords, os.path.join(self.base_name, 'data'), path)

        process_results = process_data(dataset_api=dataset_api, dataset=input_dataset,
            return_types=['embedding'], embedding_dimensions=dimensions, embedding_measures=measures,
            basis=basis, nbins=nbins, agg_function=summary)
        result = process_results['embedding']

        # write bins and coordinates
        bin_dict = dict(index=result['bins'])
        for column in result['coordinates']:
            bin_dict[column] = result['coordinates'][column]
        write_table(bin_dict, output_directory, 'index')

        for column in result['values']:
            write_table(dict(value=result['values'][column]),
                output_directory, column)

        for i in range(0, self.adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            logger.info('{} embedding {}-{}'.format(basis_name, i, end))
            process_results = process_data(dataset_api=dataset_api, dataset=input_dataset,
                return_types=['embedding'], embedding_dimensions=[], embedding_measures=adata.var_names[i:end],
                basis=basis, nbins=nbins, agg_function=summary)
            result = process_results['embedding']
            for column in result['values']:
                write_table(dict(value=result['values'][column]),
                    output_directory, column)

    def schema(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        output_file = os.path.join(self.base_name, 'index.json')
        result = SimpleData.schema(self.adata)
        del result['embeddings']
        result['summary'] = {'nObs': self.adata.shape[0], 'embeddings': []}
        for basis_name in basis_list:
            result['summary']['embeddings'].append({'basis': basis_name, 'nbins': nbins, 'agg': bin_agg_function})
        with open(output_file, 'wt') as f:
            json.dump(result, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prepare a dataset by binning on a grid using an embedding, generating feature statistics feature statistics within a category, and saving the data for easy slicing by feature.')
    parser.add_argument('dataset', help='Path to a h5ad file')

    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--basis', help='List of embeddings to precompute', action='append')
    parser.add_argument('--nbins', help='Number of bins', default=500, type=int)
    parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    args = parser.parse_args()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    logger.info('preparing ' + args.dataset + '...')
    prepare_data = PrepareData(args.dataset, args.backed, args.basis, args.nbins, args.summary)
    prepare_data.execute()

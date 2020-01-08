import argparse
import json
import logging
import os

import anndata
import cirro.data_processing as data_processing
import pandas as pd
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import EmbeddingAggregator, get_basis
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset
from cirro.io import write_table, save_adata
from cirro.simple_data import SimpleData

logger = logging.getLogger("cirro")


class PrepareData:

    def __init__(self, input_path, backed, basis_list, nbins, bin_agg_function):
        self.input_path = input_path
        self.adata = anndata.read(input_path, backed=backed)
        if self.adata.shape[0] < 50000:
            self.column_batch_size = 10000
        elif self.adata.shape[0] < 150000:
            self.column_batch_size = 5000
        else:
            self.column_batch_size = 1000

        logger.info('{} - {} x {}'.format(input_path, self.adata.shape[0], self.adata.shape[1]))
        dimensions = []
        measures = []
        others = []
        if basis_list is None:
            basis_list = list(self.adata.obsm_keys())
            # for key in basis_list:
            #     if self.adata.obsm[key].shape[1] >= 10:
            #         logger.info('Skipping {} for embedding visualization'.format(key))
            #         basis_list.remove(key)
            logger.info('Basis {}'.format(basis_list))
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
        for basis_name in basis_list:
            self.grid_embedding(basis_name, bin_agg_function, nbins)
        save_adata(self.adata, os.path.join(self.base_name, 'data'), column_batch_size=self.column_batch_size)
        self.summary_stats()
        self.grouped_stats()

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
        process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
            measures=measures, dimensions=dimensions)

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
            process_results = data_processing.handle_stats(dataset_api=dataset_api, dataset=input_dataset,
                measures=adata.var_names[i:end], dimensions=[])
            summary = process_results['summary']
            for column in summary:
                measure_summary = summary[column]
                write_table(dict(min=[measure_summary['min']],
                    max=[measure_summary['max']],
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
        logger.info('Grouped stats for {}'.format(dimensions))
        for i in range(0, self.adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            logger.info('Grouped stats {}-{}'.format(i, end))
            process_results = data_processing.handle_grouped_stats(dataset_api=dataset_api, dataset=input_dataset,
                measures=adata.var_names[i:end], dimensions=dimensions)

            dotplot_results = process_results['dotplot']

            for dotplot_result in dotplot_results:
                output_directory = os.path.join(self.base_name, 'grouped_statistics', dotplot_result['name'])
                write_table(dict(index=dotplot_result['categories']), output_directory, 'index')
                values = dotplot_result['values']
                for value in values:
                    write_table(dict(mean=value['mean'],
                        fractionExpressed=value['fractionExpressed']),
                        output_directory, value['name'])

    def grid_embedding(self, basis_name, summary, nbins):
        logger.info('{} embedding'.format(basis_name))
        basis = get_basis(basis_name, nbins, summary, False)

        dimensions = self.dimensions
        measures = self.measures
        dataset_api = self.dataset_api
        input_dataset = self.input_dataset
        adata = self.adata
        column_batch_size = self.column_batch_size
        full_basis_name = basis['full_name']

        # write cell level bin to /data
        df_with_coords = pd.DataFrame()
        obsm = adata.obsm[basis_name]
        for i in range(len(basis['coordinate_columns'])):
            df_with_coords[basis['coordinate_columns'][i]] = obsm[:, i]
        EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'],
            bin_name=full_basis_name)
        dict_with_coords = df_with_coords.to_dict(orient='list')
        write_table(dict_with_coords, os.path.join(self.base_name, 'data'), full_basis_name)
        if nbins <= 0:
            return
        # write bins and coordinates to obsm_summary/name
        result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
            measures=measures + ['__count'], dimensions=dimensions)
        output_directory = os.path.join(self.base_name, 'obsm_summary', full_basis_name)
        bin_dict = dict(index=result['bins'])
        for column in result['coordinates']:
            bin_dict[column] = result['coordinates'][column]
        write_table(bin_dict, output_directory, 'index')

        for column in result['values']:
            vals = result['values'][column]
            if not isinstance(vals, dict):
                vals = dict(value=vals)
            write_table(vals, output_directory, column)
        logger.info('{} embedding finished writing obs'.format(basis_name))
        # write X to obsm_summary/name
        for i in range(0, self.adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            logger.info('{} embedding X {}-{}'.format(basis_name, i, end))
            result = data_processing.handle_embedding(dataset_api=dataset_api, dataset=input_dataset, basis=basis,
                measures=adata.var_names[i:end], dimensions=[])
            for column in result['values']:
                write_table(dict(value=result['values'][column]),
                    output_directory, column)

    def schema(self):
        basis_list = self.basis_list
        nbins = self.nbins
        bin_agg_function = self.bin_agg_function
        output_file = os.path.join(self.base_name, 'index.json')
        result = SimpleData.schema(self.adata)
        result['precomputed'] = True
        if nbins > 0:
            result['embeddings'] = []
            for basis_name in basis_list:
                result['embeddings'].append(
                    {'precomputed': True, 'name': basis_name, 'nbins': nbins, 'agg': bin_agg_function, 'dimensions': 2})
        with open(output_file, 'wt') as f:
            json.dump(result, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prepare a dataset by binning on a grid using an embedding, generating feature statistics feature statistics within a category, and saving the data for easy slicing by feature.')
    parser.add_argument('dataset', help='Path to a h5ad file')

    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--basis', help='List of embeddings to precompute', action='append')
    parser.add_argument('--nbins', help='Number of bins. Set to 0 to disable binning', default=500, type=int)
    parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    args = parser.parse_args()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    logger.info('preparing ' + args.dataset + '...')
    prepare_data = PrepareData(args.dataset, args.backed, args.basis, args.nbins, args.summary)
    prepare_data.execute()

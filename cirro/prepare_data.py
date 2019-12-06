import argparse
import json
import logging
import os

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse
from natsort import natsorted

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis, EmbeddingAggregator
from cirro.entity import Entity
from cirro.h5ad_dataset import H5ADDataset
from cirro.parquet_dataset import ParquetDataset
from cirro.simple_data import SimpleData

logger = logging.getLogger("cirro")


def write_table(d, output_dir, name, write_statistics=True, row_group_size=None):
    os.makedirs(output_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pydict(d), output_dir + os.path.sep + name + '.parquet',
        write_statistics=write_statistics, row_group_size=row_group_size)


def fraction_expressed(x):
    return (x > 0).sum() / len(x)


class PrepareData:

    def __init__(self, input_path, backed):
        self.input_path = input_path
        self.adata = anndata.read(input_path, backed='r' if backed else None)
        dimensions = []
        measures = []
        others = []
        for i in range(len(self.adata.obs.columns)):
            name = self.adata.obs.columns[i]
            c = self.adata.obs[name]
            if pd.api.types.is_categorical_dtype(c):
                if 1 < len(c.cat.categories) < 5000:
                    dimensions.append(name)
                else:
                    others.append(name)
            elif not pd.api.types.is_string_dtype(c) and not pd.api.types.is_object_dtype(c):
                measures.append(name)
            else:
                others.append(name)
        self.others = others
        self.dimensions = dimensions
        self.measures = measures
        self.base_name = os.path.splitext(os.path.basename(input_path))[0]

    def summary_stats(self, column_batch_size=500):
        logger.info('writing summary statistics')
        # compute min, max, mean, sum for X, data.obs
        # compute value counts for data.obs categoricals
        adata = self.adata
        dimensions = self.dimensions
        measures = self.measures
        output_directory = os.path.join(self.base_name, 'counts')

        for i in range(len(dimensions)):
            name = dimensions[i]
            result = adata.obs[name].value_counts()
            logging.info('cat obs summary: {}/{}'.format(i + 1, len(dimensions)))
            write_table(dict(index=result.index, value=result.values),
                output_directory, name, write_statistics=False)
        # measure_summary_df has stats on index, measures on columns
        output_directory = os.path.join(self.base_name, 'statistics')
        measure_summary_df = adata.obs[measures].agg(['min', 'max', 'sum', 'mean'])
        for c in measure_summary_df.columns:
            write_table(dict(min=[measure_summary_df.loc['min', c]], max=[measure_summary_df.loc['max', c]],
                sum=[measure_summary_df.loc['sum', c]], mean=[measure_summary_df.loc['mean', c]]),
                output_directory, c, write_statistics=False)

        for i in range(0, adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            names = adata.var_names[slice(i, end)]
            X = adata.X[:, slice(i, end)]
            min_values = X.min(axis=0)
            max_values = X.max(axis=0)
            sums = X.sum(axis=0)
            means = X.mean(axis=0)
            num_expressed = (X > 0).sum(axis=0)
            if scipy.sparse.issparse(X):
                min_values = min_values.toarray().flatten()
                max_values = max_values.toarray().flatten()
                sums = sums.A1
                means = means.A1
                num_expressed = num_expressed.A1

            logging.info('X summary: {}/{}'.format(i + 1, adata.shape[1]))
            for j in range(len(names)):
                write_table(
                    dict(min=min_values[[j]], max=max_values[[j]], sum=sums[[j]], mean=means[[j]],
                        num_expressed=num_expressed[[j]]),
                    output_directory, names[j], write_statistics=False)


    # stats, grouped by categories in adata.obs for dot plot
    def grouped_stats(self, column_batch_size=500):
        dimensions = self.dimensions
        n_features = self.adata.shape[1]
        adata = self.adata

        for dimension in dimensions:
            logger.info('computing grouped stats for {}'.format(dimension))
            output_directory = os.path.join(self.base_name, 'grouped_statistics', dimension)

            for i in range(0, n_features, column_batch_size):
                end = i + column_batch_size
                end = min(end, n_features)
                logger.info('grouped stats for {} {}-{}/{}'.format(dimension, i, end, n_features))
                X = adata.X[:, slice(i, end)]
                if scipy.sparse.issparse(X):
                    X = X.toarray()
                df = pd.DataFrame(index=adata.obs[dimension], data=X, columns=adata.var_names[slice(i, end)].tolist())
                df = df.groupby(df.index).agg(['mean', fraction_expressed])
                summary_df = df
                sorted_categories = natsorted(summary_df.index)
                summary_df = summary_df.loc[sorted_categories]
                write_table(dict(index=summary_df.index),
                    output_directory, 'index')
                for feature in summary_df.columns.get_level_values(0):
                    write_table(dict(mean=summary_df[(feature, 'mean')],
                        fraction_expressed=summary_df[(feature, 'fraction_expressed')]),
                        output_directory, feature)

    def save_adata_X_chunk(self, start, end):
        adata = self.adata
        X_slice = adata.X[:, slice(start, end)]
        output_directory = os.path.join(self.base_name, 'data')
        names = adata.var_names[slice(start, end)]
        logger.info('writing adata X {}-{}'.format(start, end))
        for j in range(len(names)):
            X = X_slice[:, j]
            if scipy.sparse.issparse(X):
                X = X.toarray()
            indices = np.where(X != 0)[0]
            values = X[indices]

            write_table(dict(index=indices, value=values), output_directory, names[j])

    def save_adata(self, column_batch_size=500):
        logger.info('writing adata')
        adata = self.adata
        for i in range(0, adata.shape[1], column_batch_size):
            end = i + column_batch_size
            end = min(end, adata.shape[1])
            self.save_adata_X_chunk(i, end)

        logger.info('writing adata obs')
        output_directory = os.path.join(self.base_name, 'data')
        for name in adata.obs:
            # TODO sort?, row group size?
            write_table(dict(value=adata.obs[name]),
                output_directory, name)
        write_table(dict(value=adata.obs.index), output_directory, 'id')
        logger.info('writing adata obsm')

        for name in adata.obsm.keys():
            m = adata.obsm[name]
            ndim = 2  # TODO 3d
            d = {}
            for i in range(ndim):
                d[str(i)] = m[:, i].astype('float32')
            write_table(d, output_directory, name)

    def write_schema(self, basis_names, summary, nbins):
        output_file = os.path.join(self.base_name, 'index.json')
        result = SimpleData.schema(self.adata)
        del result['embeddings']
        result['summary'] = {'nObs': self.adata.shape[0], 'embeddings': []}
        for basis_name in basis_names:
            result['summary']['embeddings'].append({'basis': basis_name, 'nbins': nbins, 'agg': summary})
        with open(output_file, 'wt') as f:
            json.dump(result, f)

    def grid_embedding(self, basis_name, summary, nbins, column_batch_size=500):
        logger.info('writing embedding for {}'.format(basis_name))
        path = basis_name + '_' + str(nbins) + '_' + str(summary)
        output_directory = os.path.join(self.base_name, 'obsm_summary', path)
        return_types = ['embedding']
        basis = get_basis(basis_name)
        dataset = Entity(self.input_path, {'name': self.base_name, 'url': self.input_path})
        dataset_api = DatasetAPI()
        dataset_api.add(['pq', 'parquet'], ParquetDataset())
        ds = H5ADDataset()
        ds.path_to_data[self.base_name] = self.adata
        dataset_api.add(['h5ad'], ds)
        # write bin and coords to data
        df_with_coords = SimpleData.to_df(self.adata, [], [], [], basis=basis)
        EmbeddingAggregator.convert_coords_to_bin(df_with_coords, nbins, basis['coordinate_columns'], None)
        dict_with_coords = df_with_coords.to_dict(orient='list')
        dict_with_coords['index'] = df_with_coords.index
        write_table(dict_with_coords, os.path.join(self.base_name, 'data'), path)
        all_names = self.dimensions + self.measures + self.adata.var_names.tolist()
        first_time = True
        for i in range(0, len(all_names), column_batch_size):
            end = i + column_batch_size
            end = min(end, len(all_names))
            logger.info('{} embedding {}-{}/{}'.format(basis_name, i, end, len(all_names)))
            columns_in_batch = all_names[slice(i, end)]
            embedding_measures = []
            embedding_dimensions = []

            for column in columns_in_batch:
                if column in self.dimensions:
                    embedding_dimensions.append(column)
                else:
                    embedding_measures.append(column)
            data_processing_result = process_data(dataset_api=dataset_api, dataset=dataset,
                return_types=return_types,
                basis=basis,
                nbins=nbins, embedding_measures=embedding_measures,
                embedding_dimensions=embedding_dimensions, agg_function=summary, dotplot_measures=[],
                dotplot_dimensions=[])

            result = data_processing_result['embedding']

            for column in result['values']:
                if column == '__count' and not first_time:
                    continue
                write_table(dict(value=result['values'][column]),
                    output_directory, column)

            # write bins and coordinates
            if first_time:
                first_time = False
                bin_dict = dict(index=result['bins'])
                for column in basis['coordinate_columns']:
                    bin_dict[column] = result['coordinates'][column]
                write_table(bin_dict, output_directory, 'index')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prepare a dataset by binning on a grid using an embedding, generating feature statistics feature statistics within a category, and saving the data for easy slicing by feature.')
    parser.add_argument('dataset', help='Path to a h5ad file')
    parser.add_argument('commands', help='List of commands to execute',
        choices=['adata', 'stats', 'grouped_stats', 'basis', 'schema'],
        action='append')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--basis', help='List of embeddings to precompute', action='append')
    parser.add_argument('--nbins', help='Number of bins', default=500, type=int)
    parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    args = parser.parse_args()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())
    commands = args.commands
    prepare_data = PrepareData(args.dataset, args.backed)
    logger.info('preparing ' + args.dataset + '...')
    if 'adata' in commands:
        prepare_data.save_adata()
    if 'stats' in commands:
        prepare_data.summary_stats()
    if 'grouped_stats' in commands:
        prepare_data.grouped_stats()
    if 'basis' in commands:
        if args.basis is None:
            obsm_keys = list(prepare_data.adata.obsm_keys())
            blacklist = ['X_pca']
            for key in obsm_keys:
                if key in blacklist or p.adata.obsm[key].shape[1] >= 10:
                    obsm_keys.remove(key)
        else:
            obsm_keys = args.basis
        for b in obsm_keys:
            prepare_data.grid_embedding(b, args.summary, args.nbins)
    if 'schema' in commands:
        prepare_data.write_schema(obsm_keys, args.summary, args.nbins)

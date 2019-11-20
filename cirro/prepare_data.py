import argparse
import math
import os
import tempfile

import pandas as pd
import pyarrow
import pyarrow as pa
import pyarrow.parquet as pq
from natsort import natsorted

from cirro.data_processing import process_data
from cirro.dataset_api import DatasetAPI
from cirro.embedding_aggregator import get_basis
from cirro.entity import Entity
from cirro.parquet_dataset import ParquetDataset


class PrepareData:

    def __init__(self, input_path):
        self.input_path = input_path
        self.parquet_file = pq.ParquetFile(input_path)
        self.schema = self.parquet_file.metadata.schema.to_arrow_schema()
        measures = []
        dimensions = []
        for i in range(len(self.schema.names)):
            name = self.schema.names[i]
            data_type = self.schema.types[i]
            if isinstance(data_type, pyarrow.lib.DictionaryType):
                dimensions.append(name)
            elif not pyarrow.types.is_string(data_type):
                measures.append(name)
        self.measures = measures
        self.dimensions = dimensions

    def summary_stats(self):
        # compute min, max, mean, sum for measures
        # compute value counts for categories
        schema = self.schema
        parquet_file = self.parquet_file
        input_path = self.input_path
        measure_df = pd.DataFrame(
            index=['min', 'max', 'sum', 'mean', 'num_expressed'])  # rows are min, max, etc, columns are measures
        for i in range(len(schema.names)):
            data_type = schema.types[i]
            name = schema.names[i]
            if isinstance(data_type, pyarrow.lib.DictionaryType):
                df = parquet_file.read([name]).to_pandas()
                result = df[name].value_counts()
                result = pd.DataFrame(index=result.index, data=dict(counts=result.values))
                table = pa.Table.from_pandas(result)
                pq.write_table(table,
                    os.path.splitext(os.path.basename(input_path))[0] + '_counts_' + name + '.parquet')
            elif not pyarrow.types.is_string(data_type):
                df = parquet_file.read([name]).to_pandas()
                measure_df[name] = (df[name].min(), df[name].max(), df[name].sum(), df[name].mean(),
                                    (df[name] > 0).sum())
            else:
                print('Skipped {}'.format(name))

        pq.write_table(pa.Table.from_pandas(measure_df),
            os.path.splitext(os.path.basename(input_path))[0] + '_measure_summary.parquet')

    def dimension_stats(self, column_batch_size=1000):
        dimensions = self.dimensions
        measures = self.measures
        parquet_file = self.parquet_file
        input_path = self.input_path

        def fraction_expressed(x):
            return (x > 0).sum() / len(x)

        print('{} dimensions'.format(len(dimensions)))
        for dimension in dimensions:
            dimension_summary = None
            for i in range(0, len(measures), column_batch_size):
                end = i + column_batch_size
                end = min(end, len(measures))
                print('{} {}-{}/'.format(dimension, i, end, len(measures)))
                df = parquet_file.read([dimension] + measures[slice(i, end)]).to_pandas()
                df = df.groupby(dimension).agg(['mean', fraction_expressed])
                dimension_summary = df if dimension_summary is None else dimension_summary.join(df)
            sorted_categories = natsorted(dimension_summary.index)
            dimension_summary = dimension_summary.loc[sorted_categories]
            pq.write_table(pa.Table.from_pandas(dimension_summary),
                os.path.splitext(os.path.basename(input_path))[0] + '_statistics_' + dimension + '.parquet')

    def grid_embedding(self, basis_name, summary, nbins, column_batch_size, n_row_groups):
        input_path = self.input_path
        full_parquet_file = self.parquet_file
        full_schema = self.schema
        output = os.path.splitext(os.path.basename(input_path))[0] + '_' + basis_name + '_' + str(nbins) + '_' + str(
            summary) + '.parquet'

        # backed = args.backed

        return_types = ['embedding']

        basis = get_basis(basis_name)
        dataset = Entity(input_path, {'name': os.path.splitext(os.path.basename(input_path))[0], 'url': input_path})
        dataset_api = DatasetAPI()
        dataset_api.add(['pq', 'parquet'], ParquetDataset())
        # dataset_api.add(['h5ad'], H5ADDataset('r' if backed else None))

        full_columns = full_schema.names
        nobs = full_parquet_file.metadata.num_rows

        if n_row_groups is None:
            n_row_groups = max(1, math.ceil(nobs / 100000))
        row_group_size = math.ceil(nobs / n_row_groups)
        n_row_groups = math.ceil(row_group_size / nobs)
        tmp_files = []
        first_batch = True

        # write parquet files in batches of column_batch_size

        for i in range(0, len(full_columns), column_batch_size):
            end = i + column_batch_size
            end = min(end, len(full_columns))
            print('{}-{}/{}'.format(i, end, len(full_columns)))
            columns_in_batch = full_columns[slice(i, end)]
            embedding_measures = []
            embedding_dimensions = []

            for column in columns_in_batch:
                data_type = full_schema.types[full_schema.get_field_index(column)]
                if isinstance(data_type, pyarrow.lib.DictionaryType):
                    embedding_dimensions.append(column)
                elif not pyarrow.types.is_string(data_type):
                    embedding_measures.append(column)
                else:
                    print('Skipped {}'.format(column))
            data_processing_result = process_data(dataset_api=dataset_api, dataset=dataset,
                return_types=return_types,
                basis=basis,
                nbins=nbins, embedding_measures=embedding_measures,
                embedding_dimensions=embedding_dimensions, agg_function=summary, dotplot_measures=[],
                dotplot_dimensions=[])

            embedding_summary = data_processing_result['embedding']
            measure_df, dimension_df = embedding_summary.collect()

            # save __count, coordinate columns if first batch, bins are measure_df.index

            for column in embedding_summary.dimensions:
                measure_df[column] = dimension_df[column].values
            columns = measure_df.columns.to_list()
            if not first_batch:
                columns.remove('__count')
                for column in basis['coordinate_columns']:
                    columns.remove(column)

            table = pa.Table.from_pandas(measure_df, columns=columns)
            first_batch = False
            _, tmp_file = tempfile.mkstemp()
            tmp_files.append(tmp_file)
            pq.write_table(table, tmp_file, row_group_size=row_group_size)

        # join files
        writer = None
        for row_group in range(n_row_groups):
            df = None
            print('Group {}/{}'.format(row_group + 1, n_row_groups))
            for tmp_file in tmp_files:
                parquet_file = pq.ParquetFile(tmp_file)
                table = parquet_file.read_row_group(row_group)
                df = table.to_pandas() if df is None else pd.concat((df, table.to_pandas()), axis=1)
            table = pa.Table.from_pandas(df)
            if writer is None:
                writer = pq.ParquetWriter(output, schema=table.schema)
            writer.write_table(table)
        writer.close()
        for tmp_file in tmp_files:
            os.unlink(tmp_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prepare a dataset by binning on a grid using an embedding, computing global feature statistics, and statistics within each category')
    parser.add_argument('dataset', help='Path to a parquet file')
    parser.add_argument('--output', help='Output file')
    # parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--basis', help='Embedding basis', action='append')
    parser.add_argument('--nbins', help='Number of bins', default=500, type=int)
    parser.add_argument('--row_groups', help='Number row groups', type=int)
    parser.add_argument('--column_batch_size', help='Column batch size', default=2000, type=int)
    parser.add_argument('--summary', help='Bin summary statistic for numeric values', default='max')
    args = parser.parse_args()

    p = PrepareData(args.dataset)
    p.dimension_stats()
    p.summary_stats()
    for b in args.basis:
        p.grid_embedding(b, args.summary, args.nbins, args.column_batch_size, args.row_groups)

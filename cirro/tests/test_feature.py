import pandas as pd
from pytest import approx

from cirro.data_processing import handle_stats


def diff_measures(test_data, summary, fields):
    my_mean = test_data[:, fields].X.mean(axis=0)
    my_min = test_data[:, fields].X.min(axis=0)
    my_max = test_data[:, fields].X.max(axis=0)
    my_sum = test_data[:, fields].X.sum(axis=0)
    my_expressed = (test_data[:, fields].X != 0).sum(axis=0)

    for i in range(len(fields)):
        field = fields[i]
        field_summary = summary[field]
        assert my_min[i] == field_summary['min'], 'min'
        assert my_max[i] == field_summary['max'], 'max'
        assert my_expressed[i] == field_summary['numExpressed'], 'numExpressed'
        assert my_sum[i] == approx(field_summary['sum'], abs=0.0001), 'sum'
        assert my_mean[i] == approx(field_summary['mean'], abs=0.0000001), 'mean'


def diff_dimensions(test_data, summary, fields):
    for key in fields:
        value_counts = test_data.obs[key].value_counts()
        dimension_summary = summary[key]
        df = pd.DataFrame.from_dict(dimension_summary)
        df['categories'] = df['categories'].astype('category')
        df = df.set_index('categories')
        df = df.iloc[df.index.get_indexer_for(value_counts.index)]
        assert (df['counts'] - value_counts).sum() == 0, key


def test_measure(dataset_api, input_dataset, test_data, measures):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset, measures=measures)
    summary = process_results['summary']
    diff_measures(test_data, summary, measures)


def test_dimension(dataset_api, input_dataset, test_data, dimensions):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset,
        dimensions=dimensions, )
    summary = process_results['summary']
    diff_dimensions(test_data, summary, dimensions)


def test_measure_and_dimension(dataset_api, input_dataset, test_data, measures, dimensions):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset, measures=measures,
        dimensions=dimensions, )
    summary = process_results['summary']
    diff_measures(test_data, summary, measures)
    diff_dimensions(test_data, summary, dimensions)

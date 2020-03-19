import pandas as pd
from pytest import approx

from cirrocumulus.data_processing import handle_stats
from cirrocumulus.entity import Entity
from cirrocumulus.feature_aggregator import FeatureAggregator
from cirrocumulus.prepare_data import make_ordered, write_X_stats, write_obs_stats


def diff_measures(test_data, summary, fields, is_x):
    X = test_data[:, fields].X if is_x else test_data.obs[fields].values
    my_mean = X.mean(axis=0)
    my_min = X.min(axis=0)
    my_max = X.max(axis=0)
    my_sum = X.sum(axis=0)
    my_expressed = (X != 0).sum(axis=0)

    for i in range(len(fields)):
        field = fields[i]
        field_summary = summary[field]
        assert my_min[i] == field_summary['min'], 'min'
        assert my_max[i] == field_summary['max'], 'max'
        if 'numExpressed' in field_summary:
            assert my_expressed[i] == field_summary['numExpressed'], 'numExpressed'
        assert my_sum[i] == approx(field_summary['sum'], abs=0.001), 'sum'
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
    diff_measures(test_data, summary, measures, True)


def test_continuous_obs(dataset_api, input_dataset, test_data, continuous_obs):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset,
        measures=list(map(lambda x: 'obs/' + x, continuous_obs)))
    summary = process_results['summary']
    diff_measures(test_data, summary, continuous_obs, False)


def test_dimension(dataset_api, input_dataset, test_data, dimensions):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset,
        dimensions=dimensions, )
    summary = process_results['summary']
    diff_dimensions(test_data, summary, dimensions)


def test_all(dataset_api, input_dataset, test_data, measures, dimensions, continuous_obs):
    process_results = handle_stats(dataset_api=dataset_api, dataset=input_dataset,
        measures=measures + list(map(lambda x: 'obs/' + x, continuous_obs)),
        dimensions=dimensions, )
    summary = process_results['summary']
    diff_measures(test_data, summary, measures, True)
    diff_measures(test_data, summary, continuous_obs, False)
    diff_dimensions(test_data, summary, dimensions)


def test_save_stats(tmp_path, dataset_api, test_data, measures, continuous_obs, dimensions):
    import zarr
    url = str(tmp_path) + '.zarr'
    make_ordered(test_data.obs, None)
    test_data = test_data[:, measures]
    test_data.write_zarr(url)
    store = zarr.open(url)

    write_X_stats(store, FeatureAggregator(var_measures=measures).execute(test_data))
    write_obs_stats(store, FeatureAggregator(dimensions=dimensions, obs_measures=continuous_obs).execute(test_data))

    meta = {'name': tmp_path, 'url': url}
    input_dataset = Entity(url, meta)
    result = dataset_api.read_precomputed_stats(input_dataset, obs_keys=dimensions + continuous_obs, var_keys=measures)
    diff_measures(test_data, result, measures, True)
    diff_measures(test_data, result, continuous_obs, False)
    diff_dimensions(test_data, result, dimensions)

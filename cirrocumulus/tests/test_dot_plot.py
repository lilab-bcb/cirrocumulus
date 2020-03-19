import numpy as np
import pandas as pd
import scipy.sparse

from cirrocumulus.data_processing import handle_grouped_stats
from cirrocumulus.dotplot_aggregator import DotPlotAggregator
from cirrocumulus.entity import Entity
from cirrocumulus.prepare_data import write_grouped_stats, make_ordered


def diff_results(summarized_df, measures, results):
    dotplot_result = results[0]
    values = dotplot_result['values']
    for key in measures:
        index = -1
        for i in range(len(values)):
            if values[i]['name'] == key:
                index = i
                break
        if index == -1:
            raise ValueError(key + ' not found')
        np.testing.assert_allclose(summarized_df[key]['mean'].values, values[index]['mean'],
            atol=0.000001, err_msg='{} mean'.format(key))
        np.testing.assert_allclose(summarized_df[key]['fraction_expressed'].values,
            values[index]['fractionExpressed'], err_msg='{} fractionExpressed'.format(key))


def get_df(test_data, measures, by):
    X = test_data[:, measures].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    df = pd.DataFrame(data=X, columns=measures)
    df[by] = test_data.obs[by].values
    make_ordered(df, None)

    def fraction_expressed(g):
        return (g != 0).sum() / len(g)

    summarized_df = df.groupby(by).agg(['mean', fraction_expressed])
    return summarized_df


def test_dot_plot(dataset_api, input_dataset, test_data, measures, by):
    process_results = handle_grouped_stats(dataset_api=dataset_api, dataset=input_dataset, measures=measures,
        dimensions=[by])
    dotplot_result = process_results['dotplot']
    diff_results(get_df(test_data, measures, by), measures, dotplot_result)


def test_save_grouped_stats(tmp_path, dataset_api, test_data, measures, by):
    import zarr
    url = str(tmp_path) + '.zarr'
    make_ordered(test_data.obs, None)
    test_data = test_data[:, measures]
    test_data.write_zarr(url)
    test_data[:, measures].write_zarr(url)
    store = zarr.open(url)
    results = DotPlotAggregator(var_measures=measures, dimensions=[by]).execute(test_data)
    write_grouped_stats(store, test_data, results)
    input_dataset = Entity(url, {'name': tmp_path, 'url': url})
    result = dataset_api.read_precomputed_grouped_stats(input_dataset, obs_keys=[by], var_keys=measures)
    diff_results(get_df(test_data, measures, by), measures, result)

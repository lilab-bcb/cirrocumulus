import numpy as np
import pandas as pd
import scipy.sparse
from natsort import natsorted

from cirro.data_processing import handle_grouped_stats


def test_dot_plot(dataset_api, input_dataset, test_data, measures, by):
    X = test_data[:, measures].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    df = pd.DataFrame(data=X, columns=measures)
    df[by] = test_data.obs[by].values

    def fraction_expressed(g):
        return (g != 0).sum() / len(g)

    summarized_df = df.groupby(by).agg(['mean', fraction_expressed])
    sorted_categories = natsorted(summarized_df.index)
    summarized_df = summarized_df.loc[sorted_categories]
    process_results = handle_grouped_stats(dataset_api=dataset_api, dataset=input_dataset, measures=measures,
        dimensions=[by])
    dotplot_result = process_results['dotplot']
    dotplot_result = dotplot_result[0]
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
            atol=0.00001, err_msg='mean')
        np.testing.assert_allclose(summarized_df[key]['fraction_expressed'].values,
            values[index]['fractionExpressed'], err_msg='fractionExpressed')

import numpy as np
import pandas as pd
from flask import Blueprint, Response, request
from natsort import natsorted

from .auth_api import AuthAPI
from .database_api import DatabaseAPI
from .dataset_api import DatasetAPI
from .invalid_usage import InvalidUsage
from .util import *

blueprint = Blueprint('blueprint', __name__)

dataset_api = DatasetAPI()

auth_api = AuthAPI()
database_api = DatabaseAPI()


def check_bin_input(nbins):
    if nbins is not None:
        nbins = int(nbins)
        nbins = min(1000, nbins)
        if nbins <= 0:
            nbins = None
    return nbins


def get_basis(basis):
    if basis is not None:
        embedding_ndim = 2
        if basis.endswith('3d'):
            basis = basis[0:len(basis) - 3]
            embedding_ndim = 3
        coordinate_columns = []
        for i in range(embedding_ndim):
            coordinate_columns.append(basis + '_' + str(i + 1))
        return {'name': basis, 'dimensions': embedding_ndim, 'coordinate_columns': coordinate_columns}


@blueprint.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    return Response(error.message, error.status_code)


@blueprint.route('/server', methods=['GET'])
def handle_server():
    # no login required
    server = database_api.server()
    server['clientId'] = auth_api.client_id
    return to_json(server)


@blueprint.route('/schema', methods=['GET'])
def handle_schema():
    """Get dataset schema.
    """
    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400
    dataset = database_api.get_dataset(email, dataset_id)
    schema = dataset_api.schema(dataset)
    embeddings = schema['embeddings']
    additional_embeddings = []
    for embedding in embeddings:
        if embedding['dimensions'] == 3:
            additional_embeddings.append({'name': embedding['name'], 'dimensions': 2})
            embedding['name'] += ' 3d'
    embeddings += additional_embeddings
    embeddings = sorted(embeddings, key=lambda x: (x['name'], x['dimensions']))
    schema['embeddings'] = embeddings
    # custom_annotations = database_api.list_annotations(dataset_id)
    # schema['custom_annotations'] = custom_annotations
    return to_json(schema)


@blueprint.route('/user', methods=['GET'])
def handle_user():
    email = auth_api.auth()['email']
    user = database_api.user(email)
    return to_json(user)


def __convert_coords(df, nbins, coordinate_columns, coordinate_column_to_range=None):
    # replace coordinates with bin
    for name in coordinate_columns:
        values = df[name].values
        view_column_range = coordinate_column_to_range.get(name,
            None) if coordinate_column_to_range is not None else None
        column_min = values.min() if view_column_range is None else view_column_range[0]
        column_max = values.max() if view_column_range is None else view_column_range[1]
        df[name] = np.floor(np.interp(values, [column_min, column_max], [0, nbins - 1])).astype(int)
    if len(coordinate_columns) == 2:
        df['__bin'] = df[coordinate_columns[0]] * nbins + df[coordinate_columns[1]]
    else:
        df['__bin'] = df[coordinate_columns[2]] + nbins * (
                df[coordinate_columns[1]] + nbins * df[coordinate_columns[0]])


def __bin_df(df, nbins, coordinate_columns, reduce_function='mean', coordinate_column_to_range=None):
    __convert_coords(df, nbins, coordinate_columns, coordinate_column_to_range)
    agg_func = {}
    for column in df:
        if column not in coordinate_columns:
            if column == 'count':
                agg_func[column] = 'sum'
            elif pd.api.types.is_numeric_dtype(df[column]):
                agg_func[column] = reduce_function
            else:  # pd.api.types.is_categorical_dtype(df[column]):
                agg_func[column] = lambda x: x.mode()[0]
        else:
            agg_func[column] = 'max'
    return df.groupby('__bin').aggregate(agg_func)


def get_selected_df(basis, nbins, dataset, selectedpoints, keys, categorical_filter, index=False):
    if categorical_filter is not None:
        for key in categorical_filter:
            if key not in keys:
                keys.append(key)
    df = dataset_api.get_df(dataset, keys, basis, index=index)
    if nbins is not None:
        __convert_coords(df, nbins, basis['coordinate_columns'])
        df = df.set_index('__bin')
    filters = []
    if selectedpoints is not None:
        if nbins is not None:
            df = df[df.index.isin(selectedpoints)]
        else:
            df = df.iloc[np.array(selectedpoints)]
    if categorical_filter is not None:
        for category in categorical_filter:
            filtered_values = categorical_filter[category]
            filters.append((df[category].isin(filtered_values)))
        df = df[np.logical_and(*filters) if len(filters) > 1 else filters[0]]
    return df


def get_selected_indices_or_bins(df, nbins):
    if nbins is not None:
        indices = df.index.unique().values
        indices.sort()
        indices = indices.tolist()
    else:
        indices = df.index.values.tolist()
    return indices


def selected_value_counts(basis, nbins, dataset, selectedpoints, keys, categorical_filter):
    df = get_selected_df(basis if nbins is not None else None, nbins, dataset, selectedpoints, keys, categorical_filter)
    indices = get_selected_indices_or_bins(df, nbins)
    result = {'count': len(df), 'summary': {}}
    if nbins is not None:
        result['bins'] = indices
    else:
        result['indices'] = indices
    for column in df:
        if column not in basis['coordinate_columns'] and column != '__bin':
            if pd.api.types.is_numeric_dtype(df[column]):
                result['summary'][column] = {'num_expressed': len(df[df[column] > 0]),
                                             'variance': float(df[column].var()),
                                             'mean': float(df[column].mean())}
            else:
                result['summary'][column] = df[column].value_counts().to_dict()
    return result


@blueprint.route('/selected_ids', methods=['POST'])
def handle_selected_ids():
    email = auth_api.auth()['email']

    content = request.get_json(force=True, cache=False)
    dataset_id = content.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    dataset = database_api.get_dataset(email, dataset_id)
    basis = get_basis(content.get('embedding', None))
    keys = content['keys']
    selectedpoints = content.get('p', None)
    categorical_filter = content.get('c', None)
    # numerical_filter = content.get('n', None)
    if basis is not None:
        nbins = check_bin_input(content.get('nbins', None))
    df = get_selected_df(basis if nbins is not None else None, nbins, dataset, selectedpoints, keys, categorical_filter,
        index=True)
    return to_json(df.index.values.tolist())


@blueprint.route('/selected_value_counts', methods=['POST'])
def handle_selected_value_counts():
    email = auth_api.auth()['email']

    content = request.get_json(force=True, cache=False)
    dataset_id = content.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    dataset = database_api.get_dataset(email, dataset_id)
    basis = get_basis(content.get('embedding', None))
    keys = content['keys']
    selectedpoints = content.get('p', None)
    categorical_filter = content.get('c', None)
    if basis is not None:
        nbins = check_bin_input(content.get('nbins', None))
    result = selected_value_counts(basis, nbins, dataset, selectedpoints, keys, categorical_filter)
    return to_json(result)


@blueprint.route('/slice', methods=['GET'])
def handle_slice():
    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')
    if dataset_id == '':
        return 'Please provide an id', 400
    dataset = database_api.get_dataset(email, dataset_id)
    basis = get_basis(request.args.get('embedding', None))
    reduce_function = request.args.get('reduce_function', 'mean')
    keys = list(request.args.getlist('key'))

    if basis is not None:
        nbins = check_bin_input(request.args.get('nbins', None))

    df = dataset_api.get_df(dataset, keys, basis if basis is not None else None)

    if basis is not None and (nbins is not None or len(keys) == 0):
        df['count'] = 1.0
    result = {}
    if basis is None:  # TODO return dot plot and embedding simultaneously
        def non_zero(g):
            return np.count_nonzero(g) / g.shape[0]

        group_by_column_names = []
        for column in df:
            if not pd.api.types.is_numeric_dtype(df[column]):
                group_by_column_names.append(column)

        if len(group_by_column_names) > 1:
            df['_g'] = df[group_by_column_names].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
            df = df.drop(group_by_column_names, axis=1)
        elif len(group_by_column_names) == 1:
            df = df.rename(columns={group_by_column_names[0]: "_g"})
        else:
            df['_g'] = 'A'

        summarized_df = df.groupby('_g').aggregate([reduce_function, non_zero])
        dotplot_result = {'index': summarized_df.index.values.tolist()}
        for column in summarized_df:
            dotplot_result[column] = summarized_df[column].values.tolist()
        result['dotplot'] = dotplot_result
    if basis is not None:

        embedding_result = {'values': {}, 'summary': {}, 'coordinates': {}}

        for column in df:
            if pd.api.types.is_categorical_dtype(df[column]):
                value_counts = df[column].value_counts()
                sorted_unique_values = natsorted(value_counts.index)
                value_counts = value_counts.loc[sorted_unique_values]
                embedding_result['summary'][column] = {'values': value_counts.index.tolist(),
                                                       'counts': value_counts.values.tolist()}
            elif column not in basis['coordinate_columns'] and column != '__bin' and column != 'count':

                embedding_result['summary'][column] = {'num_expressed': len(df[df[column] > 0]),
                                                       'variance': float(df[column].var()),
                                                       'mean': float(df[column].mean())}

        if nbins is not None:
            df = __bin_df(df, nbins, basis['coordinate_columns'], reduce_function)
            # return bins to client so that client can keep track of selected bins
            embedding_result['bins'] = df.index.values.tolist()

        return_count = len(keys) == 0
        for column in df:
            if column in basis['coordinate_columns']:
                embedding_result['coordinates'][column] = df[column].values.tolist()
            elif column != '__bin':
                if column == 'count' and not return_count:
                    continue
                embedding_result['values'][column] = df[column].values.tolist()

        result['embedding'] = embedding_result

    return to_json(result)


# List available datasets
@blueprint.route('/datasets', methods=['GET'])
def handle_datasets():
    email = auth_api.auth()['email']
    return to_json(database_api.datasets(email))


@blueprint.route('/dataset', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset():
    email = auth_api.auth()['email']
    # POST=new dataset, PUT=update dataset, DELETE=delete, GET=get dataset info
    if request.method == 'PUT' or request.method == 'POST':
        content = request.get_json(force=True, cache=False)
        dataset_id = content.get('id', '')
        if request.method == 'POST' and dataset_id == '':
            return 'Please supply an id', 400
        readers = set(content.get('readers', []))
        dataset_name = content.get('name', '')
        url = content.get('url', '')  # e.g. gs://foo/a/b/c.parquet
        if dataset_name == '' or url == '':
            return 'Must supply dataset name and URL', 400
        dataset_id = database_api.upsert_dataset(email=email,
            dataset_id=dataset_id if request.method == 'PUT' else None,
            dataset_name=dataset_name, url=url, readers=readers)
        return to_json({'id': dataset_id})
    elif request.method == 'DELETE':
        content = request.get_json(force=True, cache=False)
        dataset_id = content.get('id', '')
        database_api.delete_dataset(email, dataset_id)
        return to_json('', 204)
    elif request.method == 'GET':
        dataset_id = request.args.get('id', '')
        return to_json(database_api.get_dataset(email, dataset_id, True))


    # if request.method == 'POST' or (request.method == 'PUT' and url != dataset_dict.get('url',
    #         '')):  # only check if can read on new dataset created
    # bucket = url[5:]  # remove gs://
    # slash = bucket.index('/')
    # gcp_object = urllib.parse.quote(bucket[slash + 1:], safe='')
    # bucket = urllib.parse.quote(bucket[0: slash], safe='')

    # head_request = requests.head(
    #     'https://www.googleapis.com/storage/v1/b/{bucket}/o/{object}'.format(bucket=bucket,
    #         object=gcp_object),
    #     headers={'Authorization': 'Bearer ' + request_util.credentials.get_access_token().access_token})
    # head_request.close()
    # if head_request.status_code != 200:
    #     return 'Not authorized to read {}'.format(url), 403

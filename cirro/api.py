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
        return {'name': basis, 'dimensions': embedding_ndim}


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
    schema = dataset_api.schema(dataset['url'])
    embeddings = schema['embeddings']
    additional_embeddings = []
    for embedding in embeddings:
        if embedding['dimensions'] == 3:
            additional_embeddings.append({'name': embedding['name'], 'dimensions': 2})
            embedding['name'] += ' 3d'
    embeddings += additional_embeddings
    embeddings = sorted(embeddings, key=lambda x: (x['name'], x['dimensions']))
    schema['embeddings'] = embeddings
    return to_json(schema)


@blueprint.route('/user', methods=['GET'])
def handle_user():
    email = auth_api.auth()['email']
    user = database_api.user(email)
    return to_json(user)


def __bin_coords(df, nbins, coordinate_columns, coordinate_column_to_range=None):
    # replace coordinates with bin
    for view_column_name in coordinate_columns:  # add view column _bin
        values = df[view_column_name].values
        view_column_range = coordinate_column_to_range.get(view_column_name,
            None) if coordinate_column_to_range is not None else None
        column_min = values.min() if view_column_range is None else view_column_range[0]
        column_max = values.max() if view_column_range is None else view_column_range[1]
        df[view_column_name] = np.floor(
            np.interp(values, [column_min, column_max], [0, nbins - 1])).astype(int)


def __bin(df, nbins, coordinate_columns, reduce_function, coordinate_column_to_range=None):
    __bin_coords(df, nbins, coordinate_columns, coordinate_column_to_range)
    agg_func = {}
    for column in df:
        if column not in coordinate_columns:
            if column == 'count':
                agg_func[column] = 'sum'
            elif pd.api.types.is_numeric_dtype(df[column]):
                agg_func[column] = reduce_function
            else:  # pd.api.types.is_categorical_dtype(df[column]):
                agg_func[column] = lambda x: x.mode()[0]
    if 'index' in df:
        agg_func['index'] = lambda x: x.values.tolist()
    return df.groupby(coordinate_columns, as_index=False).agg(agg_func), df[coordinate_columns]


def selected_value_counts(basis, nbins, url, selectedpoints, keys, categorical_filter):
    df = dataset_api.get_df(url, keys, basis, binary=True)
    coordinate_columns = []
    for i in range(basis['dimensions']):
        coordinate_columns.append(basis['name'] + '_' + str(i + 1))

    if nbins is not None:

        bin_df = df[coordinate_columns]
        # get indices of selected bins
        __bin_coords(bin_df, nbins, coordinate_columns, None)
        g = bin_df.groupby(coordinate_columns)
        coord_keys = g.groups.keys()
        if selectedpoints is not None and len(selectedpoints) > 0:
            coord_keys = coord_keys[selectedpoints]
        indices = []
        for coord_key in coord_keys:
            indices.append(g[coord_key].values)
        df = df.iloc[np.hstack(indices)]
    else:
        if selectedpoints is not None and len(selectedpoints) > 0:
            df = df.iloc[selectedpoints]
    if categorical_filter is not None:
        for category in categorical_filter:
            filtered_values = categorical_filter[category]
            df = df[~(df[category].isin(filtered_values))]
    result = {'count': len(df), 'indices': df.index.values.tolist(), 'categories': {}}

    for column in df:
        if column not in coordinate_columns:
            value_counts = df[column].value_counts()
            # value_count_dict = {}
            # for i in range(len(value_counts)):
            #     value_count_dict[value_counts.index[i]] = value_counts.values[i]
            result['categories'][column] = value_counts.to_dict()
    return result


# @blueprint.route('/selected_points', methods=['GET'])
# def handle_selected_points():
#     email = auth_api.auth()['email']
#     dataset_id = request.args.get('id', '')
#     if dataset_id == '':
#         return 'Please provide an id', 400
#
#     dataset = database_api.get_dataset(email, dataset_id)
#     url = dataset['url']
#     basis = get_basis(request.args.get('embedding', None))
#     values = [request.args.get('v')]
#     nbins = check_bin_input(request.args.get('nbins', None))
#     key = request.args.get('key')
#     df = dataset_api.get_df(url, [key], basis)
#
#     coordinate_columns = []
#     for i in range(basis['dimensions']):
#         coordinate_columns.append(basis['name'] + '_' + str(i + 1))
#
#     if nbins is not None:
#
#         # get indices of selected bins
#         __bin_coords(df, nbins, coordinate_columns, None)
#         g = df.groupby(coordinate_columns)
#         keys = g.groups.keys()
#         indices = []
#         for key in keys:
#             if key in values:
#                 indices.append(g[key].values)
#         result = np.hstack(indices)
#         return to_json(result.tolist())
#     else:
#         df = df.reset_index()
#         result = df[df[key].isin(values)].index.values
#         return to_json(result.index.values.tolist())


@blueprint.route('/selected_value_counts', methods=['POST'])
def handle_selected_value_counts():
    email = auth_api.auth()['email']

    content = request.get_json(force=True, cache=False)
    dataset_id = content.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    dataset = database_api.get_dataset(email, dataset_id)
    url = dataset['url']
    basis = get_basis(content.get('embedding', None))
    keys = content['keys']
    selectedpoints = content.get('p', None)
    categorical_filter = content.get('c', None)
    if basis is not None:
        nbins = check_bin_input(content.get('nbins', None))
    result = selected_value_counts(basis, nbins, url, selectedpoints, keys, categorical_filter)
    return to_json(result)


# @blueprint.route('/annotate', methods=['POST'])
# def handle_annotation():
#     email = auth_api.auth()['email']
#     content = request.get_json(force=True, cache=False)
#     dataset_id = content.get('id', '')
#     if dataset_id == '':
#         return 'Please provide an id', 400
#     dataset = database_api.get_dataset(email, dataset_id)
#     url = dataset['url']
#     ids = get_selected_cell_ids(embedding=content['embedding'], nbins=content.get('nbins', None),
#         selectedpoints=content['selectedpoints'], url=url)
#     return to_json({'id': dataset_id}), 200, get_json_headers()


@blueprint.route('/slice', methods=['GET'])
def handle_slice():
    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')
    if dataset_id == '':
        return 'Please provide an id', 400

    dataset = database_api.get_dataset(email, dataset_id)
    url = dataset['url']
    basis = get_basis(request.args.get('embedding', None))
    reduce_function = request.args.get('reduce_function', 'mean')
    keys = list(request.args.getlist('key'))

    if basis is not None:
        nbins = check_bin_input(request.args.get('nbins', None))

    df = dataset_api.get_df(url, keys, basis if basis is not None else None)

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
        return_count = len(keys) == 0
        coordinate_columns = []
        for i in range(basis['dimensions']):
            coordinate_columns.append(basis['name'] + '_' + str(i + 1))
        embedding_result = {'values': {}, 'categories': {}, 'obs': {}}
        for column in df:
            if pd.api.types.is_categorical_dtype(df[column]):
                value_counts = df[column].value_counts()
                sorted_unique_values = natsorted(value_counts.index)
                value_counts = value_counts.loc[sorted_unique_values]
                embedding_result['categories'][column] = {'values': value_counts.index.tolist(),
                                                          'counts': value_counts.values.tolist()}
            elif column not in coordinate_columns:
                value_counts = (df[column] > 0).astype('category').value_counts()  # TODO custom cut point
                embedding_result['obs'][column] = {'values': value_counts.index.tolist(),
                                                   'counts': value_counts.values.tolist()}
        if nbins is not None:
            df, df_with_coords = __bin(df, nbins, coordinate_columns, reduce_function)

        for column in df:
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

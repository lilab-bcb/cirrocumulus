from flask import Blueprint, Response, request

import cirrocumulus.data_processing as data_processing
from cirrocumulus.embedding_aggregator import get_basis
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


@blueprint.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    return Response(error.message, error.status_code)


@blueprint.route('/server', methods=['GET'])
def handle_server():
    # no login required
    server = database_api.server()
    server['clientId'] = auth_api.client_id
    return to_json(server)


@blueprint.route('/filters', methods=['GET'])
def handle_dataset_filters():
    """List filters available for a dataset.
    """
    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    database_api.get_dataset(email, dataset_id)
    dataset_filters = database_api.dataset_filters(email, dataset_id)
    # {'id': result.id, 'name': result['name']}
    return to_json(dataset_filters)


@blueprint.route('/export_filters', methods=['GET'])
def handle_export_dataset_filters():
    """Download filters in a csv file for a dataset.
    """
    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    dataset = database_api.get_dataset(email, dataset_id)
    dataset_filters = database_api.dataset_filters(email, dataset_id)
    text = data_processing.handle_export_dataset_filters(dataset_api, dataset, dataset_filters)
    return Response(text, mimetype='text/plain')


@blueprint.route('/filter', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset_filter():
    """CRUD for a dataset filter.
    """

    email = auth_api.auth()['email']
    if request.method == 'GET':
        dataset_filter_id = request.args.get('id', '')
        if dataset_filter_id == '':
            return 'Please provide an id', 400
        return to_json(database_api.get_dataset_filter(email, dataset_filter_id))
    content = request.get_json(force=True, cache=False)
    filter_id = content.get('id')
    dataset_id = content.get('ds_id')
    # POST=new, PUT=update , DELETE=delete, GET=get
    if request.method == 'PUT' or request.method == 'POST':
        if request.method == 'PUT' and filter_id is None:
            return 'Please supply an id', 400
        if request.method == 'POST' and dataset_id is None:
            return 'Please supply a ds_id', 400
        dataset_filter = content.get('value')
        filter_name = content.get('name')
        filter_notes = content.get('notes')
        filter_id = database_api.upsert_dataset_filter(
            email=email,
            dataset_id=dataset_id,
            filter_id=filter_id if request.method == 'PUT' else None,
            dataset_filter=dataset_filter,
            filter_name=filter_name,
            filter_notes=filter_notes)
        return to_json({'id': filter_id})
    elif request.method == 'DELETE':
        database_api.delete_dataset_filter(email, filter_id)
        return to_json('', 204)


def get_schema_and_dataset():
    """Get dataset schema and dataset object
    """

    email = auth_api.auth()['email']
    dataset_id = request.args.get('id', '')
    if dataset_id == '':
        return 'Please provide an id', 400
    dataset = database_api.get_dataset(email, dataset_id)
    schema = dataset_api.schema(dataset)
    return schema, dataset


@blueprint.route('/schema', methods=['GET'])
def handle_schema():
    """Get dataset schema.
    """
    return to_json(get_schema_and_dataset()[0])


# @blueprint.route('/download', methods=['GET'])
# def handle_file():
#     email = auth_api.auth()['email']
#     dataset_id = request.args.get('id', '')
#     dataset = database_api.get_dataset(email, dataset_id, False)
#     file = request.args.get('file')
#     path = dataset['url']
#     import os
#     path = os.path.dirname(path)
#
#     file_path = os.path.join(path, file)
#     chunk_size = 4096
#     with dataset_api.fs.open(file_path) as s:
#         def generate():
#             while True:
#                 chunk = s.read(chunk_size)
#                 if len(chunk) <= 0:
#                     break
#                 yield chunk
#
#         return Response(generate())  # TODO, correct mimetype


@blueprint.route('/user', methods=['GET'])
def handle_user():
    email = auth_api.auth()['email']
    user = database_api.user(email)
    return to_json(user)


def get_email_and_dataset(content):
    email = auth_api.auth()['email']
    dataset_id = content.get('id', '')
    if dataset_id == '':
        return 'Please supply an id', 400
    dataset = database_api.get_dataset(email, dataset_id)
    return email, dataset


# @blueprint.route('/diff_exp', methods=['POST'])
# def handle_diff_exp():
#     content = request.get_json(cache=False)
#     email, dataset = get_email_and_dataset(content)
#     data_filter = content.get('filter')
#     var_range = content.get('var_range')
#     return to_json(
#         data_processing.handle_diff_exp(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter,
#             var_range=var_range))


@blueprint.route('/embedding', methods=['POST'])
def handle_embedding():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)

    nbins = check_bin_input(content.get('nbins', None))
    agg_function = content.get('agg', 'max')
    ndim = content.get('ndim', '2')
    precomputed = content.get('precomputed', False)
    dimensions = content.get('dimensions', [])
    measures = content.get('measures', [])
    basis = get_basis(content.get('basis'), nbins=nbins, agg=agg_function, dimensions=ndim, precomputed=precomputed)
    return to_json(
        data_processing.handle_embedding(dataset_api=dataset_api, dataset=dataset, basis=basis, measures=measures,
            dimensions=dimensions))


@blueprint.route('/selection', methods=['POST'])
def handle_selection():
    # selection includes stats, coordinates
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    data_filter = content.get('filter')
    return to_json(data_processing.handle_selection(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter,
        measures=content.get('measures', []), dimensions=content.get('dimensions', []),
        embeddings=content.get('embeddings', [])))


@blueprint.route('/stats', methods=['POST'])
def handle_stats():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    measures = content.get('measures', [])
    dimensions = content.get('dimensions', [])
    return to_json(data_processing.handle_stats(dataset_api=dataset_api, dataset=dataset,
        measures=measures, dimensions=dimensions))


@blueprint.route('/grouped_stats', methods=['POST'])
def handle_grouped_stats():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    dimensions = content.get('dimensions', [])
    measures = content.get('measures', [])
    return to_json(data_processing.handle_grouped_stats(dataset_api=dataset_api, dataset=dataset,
        measures=measures, dimensions=dimensions))


@blueprint.route('/selected_ids', methods=['POST'])
def handle_selected_ids():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    data_filter = content.get('filter')
    return to_json(
        data_processing.handle_selection_ids(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter))


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

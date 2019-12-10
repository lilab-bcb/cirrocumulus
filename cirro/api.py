from cirro.data_processing import process_data
from cirro.embedding_aggregator import get_basis
from flask import Blueprint, Response, request

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

    dataset = database_api.get_dataset(email, dataset_id)
    dataset_filters = database_api.dataset_filters(email, dataset_id)
    return to_json(dataset_filters)


@blueprint.route('/filter', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset_filter():
    """CRUD for a dataset filter.
    """

    email = auth_api.auth()['email']
    if request.method == 'GET':
        dataset_filter_id = request.args.get('id', '')
        if dataset_filter_id == '':
            return 'Please provide an id', 400
        return to_json( database_api.get_dataset_filter(email, dataset_filter_id))
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
    embeddings = schema.get('embeddings')
    if embeddings is not None:
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


def _handle_slice(content):
    email = auth_api.auth()['email']
    dataset_id = content.get('id', '')
    if dataset_id is '':
        return 'Please supply an id', 400
    dataset = database_api.get_dataset(email, dataset_id)
    basis = get_basis(content.get('embedding', None))
    nbins = None
    if basis is not None:
        nbins = check_bin_input(content.get('nbins', None))
    return_types = set(content.get('types'))
    agg_function = content.get('agg', 'mean')
    embedding_dimensions = content.get('embedding_dimensions', [])
    embedding_measures = content.get('embedding_measures', [])
    dotplot_dimensions = content.get('dotplot_dimensions', [])
    dotplot_measures = content.get('dotplot_measures', [])
    summary_measures = content.get('summary_measures', [])
    summary_dimensions = content.get('summary_dimensions', [])
    data_filter = content.get('filter', None)
    return process_data(dataset_api=dataset_api, dataset=dataset, return_types=return_types, basis=basis, nbins=nbins,
        embedding_measures=embedding_measures, embedding_dimensions=embedding_dimensions,
        dotplot_measures=dotplot_measures, dotplot_dimensions=dotplot_dimensions, summary_measures=summary_measures,
        summary_dimensions=summary_dimensions, agg_function=agg_function,
        data_filter=data_filter)


@blueprint.route('/slice', methods=['POST'])
def handle_slice():
    data_processing_result = _handle_slice(request.get_json(cache=False))
    return to_json(data_processing_result)


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

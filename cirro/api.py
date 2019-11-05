import pandas as pd
from cirro.data_processing import process_data
from cirro.embedding_aggregator import get_basis
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


def _handle_slice(content):
    email = auth_api.auth()['email']
    dataset_id = content.get('id', '')
    dataset = database_api.get_dataset(email, dataset_id)
    basis = get_basis(content.get('embedding', None))
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
    return process_data(dataset_api=dataset_api, dataset=dataset, basis=basis, nbins=nbins,
        embedding_measures=embedding_measures, embedding_dimensions=embedding_dimensions,
        dotplot_measures=dotplot_measures, dotplot_dimensions=dotplot_dimensions, return_types=return_types,
        agg_function=agg_function, summary_dimensions=summary_dimensions, summary_measures=summary_measures,
        data_filter=data_filter)


@blueprint.route('/slice', methods=['POST'])
def handle_slice():
    data_processing_result = _handle_slice(request.get_json(cache=False))
    json_result = {}
    if 'summary' in data_processing_result:
        measure_feature_summary, dimensions_feature_summary = data_processing_result['summary'].collect()
        json_result['summary'] = {'measures': {}, 'dimensions': {}}
        if measure_feature_summary is not None:
            unique_features = measure_feature_summary.columns.unique(0)
            features = measure_feature_summary.columns.get_level_values(0)
            for feature in unique_features:
                feature_stat = {}
                feature_df = measure_feature_summary.iloc[:, features == feature]
                for column in feature_df:
                    statistic_name = column[1]
                    feature_stat[statistic_name] = feature_df[column].values.tolist()
                json_result['summary']['measures'][feature] = feature_stat

        for dimension_name in dimensions_feature_summary:
            dimension_df = dimensions_feature_summary[dimension_name]
            dimension_df.index = pd.MultiIndex.from_tuples(dimension_df.index)
            dimension_df = dimension_df.reset_index()
            dimension_df = dimension_df.pivot(index='level_1', columns='level_0')
            sorted_categories = natsorted(dimension_df.index)
            dimension_df = dimension_df.loc[sorted_categories]
            key = (dimension_name, True)
            dimension_summary = {'categories': dimension_df.index.values.tolist()}
            if key in dimension_df:
                dimension_summary['selected_counts'] = dimension_df[key].values.tolist()
            key = (dimension_name, False)
            if key in dimension_df:
                dimension_summary['unselected_counts'] = dimension_df[key].values.tolist()
            json_result['summary']['dimensions'][dimension_name] = dimension_summary

    if 'embedding' in data_processing_result:
        embedding_summary = data_processing_result['embedding']
        measure_df, dimension_df = embedding_summary.collect()
        embedding_json = {'coordinates': {}, 'values': {}}
        json_result['embedding'] = embedding_json
        if embedding_summary.nbins is not None:
            embedding_json['bins'] = measure_df.index.values.tolist()
        for column in embedding_summary.basis['coordinate_columns']:
            embedding_json['coordinates'][column] = measure_df[column].values.tolist()
        for column in embedding_summary.measures:
            if column == '__count' and not embedding_summary.count:
                continue
            embedding_json['values'][column] = measure_df[column].values.tolist()
        for column in embedding_summary.dimensions:
            embedding_json['values'][column] = dimension_df[
                column].values.tolist() if embedding_summary.nbins is not None else measure_df[column].values.tolist()

    if 'dotplot' in data_processing_result:
        json_result['dotplot'] = data_processing_result['dotplot'].collect()
    if 'selection' in data_processing_result:
        indices_or_bins, count = data_processing_result['selection'].collect()
        json_result['selection'] = {'indices_or_bins': indices_or_bins.values.tolist(), 'count': count}
    if 'ids' in data_processing_result:
        json_result['ids'] = data_processing_result['ids'].collect()
    return to_json(json_result)


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

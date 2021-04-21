import os

from flask import Blueprint, Response, request, stream_with_context, current_app

import cirrocumulus.data_processing as data_processing
from .dataset_api import DatasetAPI
from .envir import CIRRO_SERVE, CIRRO_FOOTER, CIRRO_UPLOAD, CIRRO_BRAND
from .file_system_adapter import get_scheme
from .invalid_usage import InvalidUsage
from .job_api import submit_job
from .util import to_json

blueprint = Blueprint('blueprint', __name__)

dataset_api = DatasetAPI()


@blueprint.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    return Response(error.message, error.status_code)


def get_database():
    return current_app.config['DATABASE']


def get_auth():
    return current_app.config['AUTH']


@blueprint.route('/server', methods=['GET'])
def handle_server():
    # no login required
    server = get_database().server()
    server['clientId'] = get_auth().client_id
    if os.environ.get(CIRRO_FOOTER) is not None:
        with open(os.environ.get(CIRRO_FOOTER), 'rt') as f:
            server['footer'] = f.read()
    if os.environ.get(CIRRO_BRAND) is not None:
        with open(os.environ.get(CIRRO_BRAND), 'rt') as f:
            server['brand'] = f.read()
    # server['brand'] = os.environ.get(CIRRO_BRAND)
    server['jobs'] = os.environ.get('GAE_APPLICATION') is None
    server['upload'] = os.environ.get(CIRRO_UPLOAD) is not None
    return to_json(server)


@blueprint.route('/filters', methods=['GET'])
def handle_dataset_filters():
    """List filters available for a dataset.
    """
    email = get_auth().auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    dataset_filters = get_database().dataset_filters(email, dataset_id)
    # {'id': result.id, 'name': result['name']}
    return to_json(dataset_filters)


@blueprint.route('/export_filters', methods=['GET'])
def handle_export_dataset_filters():
    """Download filters in a csv file for a dataset.
    """
    email = get_auth().auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400
    database_api = get_database()
    dataset_filters = database_api.dataset_filters(email, dataset_id)
    dataset = database_api.get_dataset(email, dataset_id)
    text = data_processing.handle_export_dataset_filters(dataset_api, dataset, dataset_filters)
    return Response(text, mimetype='text/plain')


@blueprint.route('/category_name', methods=['GET', 'PUT'])
def handle_category_name():
    """CRUD for category name.
    """

    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'GET':
        dataset_id = request.args.get('id', '')
        if dataset_id == '':
            return 'Please provide an id', 400
        return to_json(database_api.category_names(email, dataset_id))
    content = request.get_json(force=True, cache=False)
    category = content.get('c')
    original_name = content.get('o')
    new_name = content.get('n')
    dataset_id = content.get('id')
    if request.method == 'PUT':
        database_api.upsert_category_name(
            email=email,
            category=category,
            dataset_id=dataset_id,
            original_name=original_name,
            new_name=new_name)
        return '', 200


@blueprint.route('/feature_set', methods=['POST', 'DELETE'])
def handle_feature_set():
    """CRUD for a feature set.
    """

    email = get_auth().auth()['email']
    database_api = get_database()
    # if request.method == 'GET':
    #     filter_id = request.args.get('id', '')
    #     dataset_id = request.args.get('ds_id', '')
    #     if filter_id == '' or dataset_id == '':
    #         return 'Please provide an id', 400
    #     return to_json(database_api.get_dataset_filter(email, dataset_id=dataset_id, filter_id=filter_id))
    content = request.get_json(force=True, cache=False)
    set_id = content.get('id')
    dataset_id = content.get('ds_id')
    # POST=new, PUT=update , DELETE=delete, GET=get
    if request.method == 'PUT' or request.method == 'POST':
        if request.method == 'PUT' and set_id is None:
            return 'Please supply an id', 400
        if request.method == 'POST' and dataset_id is None:
            return 'Please supply a ds_id', 400

        name = content.get('name')
        category = content.get('category')
        features = content.get('features')
        set_id = database_api.upsert_feature_set(
            email=email,
            dataset_id=dataset_id,
            set_id=set_id if request.method == 'PUT' else None,
            name=name,
            category=category,
            features=features)
        return to_json({'id': set_id})
    elif request.method == 'DELETE':
        database_api.delete_feature_set(email, dataset_id=dataset_id, set_id=set_id)
        return to_json('', 204)


@blueprint.route('/filter', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset_filter():
    """CRUD for a dataset filter.
    """

    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'GET':
        filter_id = request.args.get('id', '')
        dataset_id = request.args.get('ds_id', '')
        if filter_id == '' or dataset_id == '':
            return 'Please provide an id', 400
        return to_json(database_api.get_dataset_filter(email, dataset_id=dataset_id, filter_id=filter_id))
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
        database_api.delete_dataset_filter(email, dataset_id=dataset_id, filter_id=filter_id)
        return to_json('', 204)


@blueprint.route('/schema', methods=['GET'])
def handle_schema():
    """Get dataset schema.
    """
    email = get_auth().auth()['email']
    database_api = get_database()
    dataset_id = request.args.get('id', '')
    dataset = database_api.get_dataset(email, dataset_id)
    schema = dataset_api.schema(dataset)
    schema.update(dataset)  # add title, etc from database to schema
    schema['markers'] = database_api.get_feature_sets(email=email, dataset_id=dataset_id)
    return to_json(schema)


def get_file_path(file, dataset_url):
    # when serving dataset, image must be relative to dataset directory
    if os.environ.get(CIRRO_SERVE) == 'true':
        _, ext = os.path.splitext(dataset_url)
        if file[0] == '/' or file.find('..') != -1:
            raise ValueError('Incorrect path')
        file_path = os.path.join(dataset_url, file)
    else:
        file_path = file
        scheme = get_scheme(file_path)
        if scheme == 'file' and not os.path.exists(file_path):
            file_path = os.path.join(dataset_url, file)
    return file_path


def send_file(file_path):
    import mimetypes
    mimetype, encoding = mimetypes.guess_type(file_path)
    # with dataset_api.fs_adapter.get_fs(file_path).open(file_path) as f:
    #     bytes = f.read()
    # return Response(bytes, mimetype=mimetype[0])
    chunk_size = 4096
    f = dataset_api.fs_adapter.get_fs(file_path).open(file_path)

    def generate():
        while True:
            chunk = f.read(chunk_size)
            if len(chunk) <= 0:
                f.close()
                break
            yield chunk

    r = Response(stream_with_context(generate()), mimetype=mimetype)
    r.headers["Content-Encoding"] = encoding
    return r


@blueprint.route('/file', methods=['GET'])
def handle_file():
    email = get_auth().auth()['email']
    database_api = get_database()
    dataset_id = request.args.get('id', '')
    dataset = database_api.get_dataset(email, dataset_id)
    file_path = get_file_path(request.args.get('file'), dataset['url'])
    return send_file(file_path)


@blueprint.route('/user', methods=['GET'])
def handle_user():
    email = get_auth().auth()['email']
    database_api = get_database()
    user = database_api.user(email)
    return to_json(user)


def get_email_and_dataset(content):
    email = get_auth().auth()['email']
    dataset_id = content.get('id', '')
    if dataset_id == '':
        return 'Please supply an id', 400
    database_api = get_database()
    dataset = database_api.get_dataset(email, dataset_id)
    return email, dataset


@blueprint.route('/data', methods=['POST'])
def handle_data():
    json_request = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(json_request)

    return to_json(
        data_processing.handle_data(dataset_api=dataset_api, dataset=dataset,
            embedding_list=json_request.get('embedding'),
            values=json_request.get('values'),
            stats=json_request.get('stats'),
            grouped_stats=json_request.get('groupedStats'),
            selection=json_request.get('selection')))


@blueprint.route('/selection', methods=['POST'])
def handle_selection():
    # selection includes stats, coordinates
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    data_filter = content.get('filter')
    return to_json(data_processing.handle_selection(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter,
        measures=content.get('measures', []), dimensions=content.get('dimensions', []),
        embeddings=content.get('embeddings', [])))


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
    email = get_auth().auth()['email']
    database_api = get_database()
    return to_json(database_api.datasets(email))


def upload_file(file):
    from werkzeug.utils import secure_filename
    filename = secure_filename(file.filename)
    dest = os.path.join(os.environ.get(CIRRO_UPLOAD), filename)
    dest_fs = dataset_api.fs_adapter.get_fs(dest).open(dest, mode='wb')
    file.save(dest_fs)
    dest_fs.close()
    return dest


def copy_url(url):
    import fsspec
    from werkzeug.utils import secure_filename
    upload = os.environ.get(CIRRO_UPLOAD)

    if url.endswith('/'):  # don't copy directories
        return url
    from urllib.parse import urlparse
    if urlparse(upload).netloc == urlparse(url).netloc:  # don't copy if already in the same bucket
        return url
    src_scheme = get_scheme(url)
    src_fs = fsspec.filesystem(src_scheme)
    filename = secure_filename(os.path.basename(url))
    dest = os.path.join(upload, filename)
    dest_scheme = get_scheme(dest)
    dest_fs = fsspec.filesystem(dest_scheme)
    out = dest_fs.open(dest, 'wb')
    n = 1024 * 1024
    with src_fs.open(url, mode='rb') as r:
        while True:
            b = r.read(n)
            if not b:
                break
            out.write(b)
    out.close()
    return dest


@blueprint.route('/dataset', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset():
    email = get_auth().auth()['email']
    database_api = get_database()
    # POST=new dataset, PUT=update dataset, DELETE=delete, GET=get dataset info
    if request.method == 'PUT' or request.method == 'POST':
        dataset_id = request.form.get('id')
        dataset_name = request.form.get('name')
        description = request.form.get('description')
        species = request.form.get('species')
        title = request.form.get('title')
        url = request.form.get('url')  # e.g. gs://foo/a/b/
        file = request.files.get('file')
        readers = request.form.get('readers')
        if readers is not None:
            import json
            readers = json.loads(readers)
        if request.method == 'PUT' and dataset_id is None:  # update
            return 'Please supply an id', 400
        if request.method == 'POST' and dataset_name == '':  # new
            return 'Must supply dataset name', 400

        if url is not None and os.environ.get(CIRRO_UPLOAD) is not None:
            url = copy_url(url)

        if file is not None:
            if os.environ.get(CIRRO_UPLOAD) is None:
                return 'Upload not supported', 400
            url = upload_file(file)
        dataset_id = database_api.upsert_dataset(email=email,
            dataset_id=dataset_id if request.method == 'PUT' else None,
            dataset_name=dataset_name, url=url, readers=readers, description=description, title=title, species=species)
        return to_json({'id': dataset_id})
    elif request.method == 'DELETE':
        content = request.get_json(force=True, cache=False)
        dataset_id = content.get('id', '')
        database_api.delete_dataset(email, dataset_id)
        return to_json('', 204)
    elif request.method == 'GET':
        dataset_id = request.args.get('id', '')
        return to_json(database_api.get_dataset(email, dataset_id, True))


@blueprint.route('/job', methods=['GET', 'DELETE'])
def handle_job():
    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'DELETE':
        content = request.get_json(force=True, cache=False)
        job_id = content.get('id', '')
        database_api.delete_job(email, job_id)
        return to_json('', 204)
    else:
        job_id = request.args.get('id', '')
        job_status = request.args.get('status', '0')
        return_result = job_status == '0'
        if job_id.startswith('cirro-'):  # precomputed result
            dataset_id = request.args.get('ds_id', '')
            email = get_auth().auth()['email']
            dataset = database_api.get_dataset(email, dataset_id)
            file_path = get_file_path(os.path.join('uns', job_id + '.json.gz'), dataset['url'])
            return send_file(file_path)
        result = database_api.get_job(email=email, job_id=job_id, return_result=return_result)
        if return_result:
            result = Response(result, mimetype='application/json')
            result.headers["Content-Encoding"] = 'gzip'
        return result


@blueprint.route('/jobs', methods=['GET'])
def handle_jobs():
    email = get_auth().auth()['email']
    database_api = get_database()
    ds_id = request.args.get('id', '')
    return to_json(database_api.get_jobs(email=email, dataset_id=ds_id))


@blueprint.route('/submit_job', methods=['POST'])
def handle_submit_job():
    if os.environ.get('GAE_APPLICATION') is None:  # TODO
        database_api = get_database()
        content = request.get_json(force=True, cache=False)
        email, dataset = get_email_and_dataset(content)
        params = content.get('params')
        job_type = content.get('type')
        job_name = content.get('name')
        return dict(id=submit_job(database_api=database_api, dataset_api=dataset_api, email=email, dataset=dataset,
            job_name=job_name, job_type=job_type, params=params))

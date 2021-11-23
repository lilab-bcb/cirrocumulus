import json
import os
from urllib.parse import urlparse

from flask import Blueprint, Response, request, stream_with_context

import cirrocumulus.data_processing as data_processing
from .anndata_util import adata_to_df
from .blueprint_util import get_database, map_url, get_auth
from .dataset_api import DatasetAPI
from .envir import CIRRO_SERVE, CIRRO_FOOTER, CIRRO_UPLOAD, CIRRO_BRAND, CIRRO_EMAIL, CIRRO_DATASET_SELECTOR_COLUMNS, \
    CIRRO_CELL_ONTOLOGY, CIRRO_STATIC_DIR, CIRRO_MIXPANEL, CIRRO_SPECIES
from .invalid_usage import InvalidUsage
from .job_api import submit_job, delete_job
from .util import json_response, get_scheme, get_fs, open_file

cirro_blueprint = Blueprint('cirro', __name__)

dataset_api = DatasetAPI()


def get_email_and_dataset(content):
    email = get_auth().auth()['email']
    dataset_id = content.get('id', '')
    if dataset_id == '':
        return 'Please supply an id', 400
    database_api = get_database()
    dataset = database_api.get_dataset(email, dataset_id)
    dataset['url'] = map_url(dataset['url'])
    return email, dataset


def upload_file(file):
    from werkzeug.utils import secure_filename
    filename = secure_filename(file.filename)
    dest = os.path.join(os.environ.get(CIRRO_UPLOAD), filename)
    dest_fs = get_fs(dest).open(dest, mode='wb')
    file.save(dest_fs)
    dest_fs.close()
    return dest


def copy_url(url):
    from werkzeug.utils import secure_filename
    upload = os.environ.get(CIRRO_UPLOAD)

    if url.endswith('/') or os.path.isdir(url):  # don't copy directories
        return url

    if urlparse(upload).netloc == urlparse(url).netloc:  # don't copy if already in the same bucket
        return url
    src_scheme = get_scheme(url)
    src_fs = get_fs(src_scheme)
    filename = secure_filename(os.path.basename(url))
    dest = os.path.join(upload, filename)
    dest_scheme = get_scheme(dest)
    dest_fs = get_fs(dest_scheme)
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


def get_file_path(file, dataset_url):
    # when serving dataset, file must be relative to dataset directory
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
    chunk_size = 4096
    f = get_fs(file_path).open(file_path)

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


@cirro_blueprint.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    return Response(error.message, error.status_code)


@cirro_blueprint.route('/server', methods=['GET'])
def handle_server():
    # no login required
    d = {}
    d['email'] = os.environ.get(CIRRO_EMAIL)
    d['capabilities'] = get_database().capabilities()
    d['clientId'] = get_auth().client_id
    if os.environ.get(CIRRO_MIXPANEL) is not None:
        d['mixpanel'] = os.environ[CIRRO_MIXPANEL]
    if os.environ.get(CIRRO_DATASET_SELECTOR_COLUMNS) is not None:
        if os.path.exists(os.environ[CIRRO_DATASET_SELECTOR_COLUMNS]):
            with open(os.environ[CIRRO_DATASET_SELECTOR_COLUMNS], 'rt') as f:
                d['datasetSelectorColumns'] = json.loads(f.read())
        else:
            d['datasetSelectorColumns'] = json.loads(os.environ[CIRRO_DATASET_SELECTOR_COLUMNS])
    if os.environ.get(CIRRO_FOOTER) is not None:
        with open(os.environ.get(CIRRO_FOOTER), 'rt') as f:
            d['footer'] = f.read()
    if os.environ.get(CIRRO_BRAND) is not None:
        with open(os.environ.get(CIRRO_BRAND), 'rt') as f:
            d['brand'] = f.read()
    d['upload'] = os.environ.get(CIRRO_UPLOAD) is not None
    if os.environ.get(CIRRO_SPECIES) is not None:
        if get_fs(os.environ[CIRRO_SPECIES]).exists(os.environ[CIRRO_SPECIES]):
            with open_file(os.environ[CIRRO_SPECIES], 'rt') as f:
                d['species'] = json.load(f)
        else:
            d['species'] = json.loads(os.environ.get(CIRRO_SPECIES))
    else:
        d['species'] = dict(favorite=["Homo sapiens", "Mus musculus"], other=["Gallus gallus", "Macaca fascicularis",
                                                                              "Macaca mulatta", "Rattus norvegicus"]);

    if os.environ.get(CIRRO_CELL_ONTOLOGY) is not None:
        if get_fs(os.environ[CIRRO_CELL_ONTOLOGY]).exists(os.environ[CIRRO_CELL_ONTOLOGY]):
            def parse_term(f):
                term = dict()
                for term_line in f:
                    term_line = term_line.strip()
                    if term_line.startswith("["):
                        break
                    index = term_line.find(':')
                    if index != -1:
                        key = term_line[:index]
                        term[key] = term_line[index + 1:].strip()
                # if term.get('is_obsolete', '') == 'true':
                #     return None
                return term

            terms = []
            with open_file(os.environ[CIRRO_CELL_ONTOLOGY], 'rt') as f:  # obo file
                for line in f:
                    line = line.strip()
                    if line == "[Term]":
                        t = parse_term(f)
                        if t is not None:
                            terms.append(t)

            d['ontology'] = dict(cellTypes=terms)
        else:
            d['ontology'] = dict(cellTypes=os.environ[CIRRO_CELL_ONTOLOGY])
    return json_response(d)


@cirro_blueprint.route('/category_name', methods=['GET', 'PUT'])
def handle_category_name():
    """CRUD for category name.
    """

    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'GET':
        dataset_id = request.args.get('id', '')
        if dataset_id == '':
            return 'Please provide an id', 400
        return json_response(database_api.category_names(email, dataset_id))

    if request.method == 'PUT':
        content = request.get_json(force=True, cache=False)
        category = content['name']
        dataset_id = content['id']
        original_value = content.get('originalValue')
        update = {}
        if 'newValue' in content:
            update['newValue'] = content['newValue']
        if 'positiveMarkers' in content:
            update['positiveMarkers'] = content['positiveMarkers']
        if 'negativeMarkers' in content:
            update['negativeMarkers'] = content['negativeMarkers']
        if 'color' in content:
            update['color'] = content['color']
        database_api.upsert_category_name(
            email=email,
            category=category,
            dataset_id=dataset_id,
            original_value=original_value,
            update=update)
        return '', 200


@cirro_blueprint.route('/feature_set', methods=['POST', 'DELETE'])
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
    #     return json_response(database_api.get_dataset_filter(email, dataset_id=dataset_id, filter_id=filter_id))
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
        return json_response({'id': set_id})
    elif request.method == 'DELETE':
        database_api.delete_feature_set(email, dataset_id=dataset_id, set_id=set_id)
        return json_response('', 204)


@cirro_blueprint.route('/views', methods=['GET'])
def handle_dataset_views():
    """List available views for a dataset.
    """
    email = get_auth().auth()['email']
    dataset_id = request.args.get('id', '')

    if dataset_id == '':
        return 'Please provide an id', 400

    return json_response(get_database().dataset_views(email, dataset_id))


@cirro_blueprint.route('/view', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset_view():
    """CRUD for a dataset view.
    """
    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'GET':
        view_id = request.args.get('id', '')
        if view_id == '':
            return 'Please provide an id', 400
        return json_response(database_api.get_dataset_view(email, view_id=view_id))
    d = request.get_json(force=True, cache=False)

    # POST=new, PUT=update , DELETE=delete, GET=get
    if request.method == 'PUT' or request.method == 'POST':
        dataset_id = d.pop('ds_id')
        view_id = d.pop('id') if request.method == 'PUT' else None
        if request.method == 'PUT' and view_id is None:
            return 'Please supply an id', 400
        if request.method == 'POST' and dataset_id is None:
            return 'Please supply a ds_id', 400

        result = database_api.upsert_dataset_view(
            email=email,
            dataset_id=dataset_id,
            view=d)
        return json_response(result)
    elif request.method == 'DELETE':
        database_api.delete_dataset_view(email, view_id=d.pop('id'))
        return json_response('', 204)


@cirro_blueprint.route('/schema', methods=['GET'])
def handle_schema():
    """Get dataset schema.
    """
    email = get_auth().auth()['email']
    database_api = get_database()
    dataset_id = request.args.get('id', '')
    dataset = database_api.get_dataset(email, dataset_id)
    dataset['url'] = map_url(dataset['url'])
    schema = dataset  # dataset has title, etc. from database
    schema['markers'] = database_api.get_feature_sets(email=email, dataset_id=dataset_id)
    schema.update(dataset_api.get_schema(dataset))
    return json_response(schema)


@cirro_blueprint.route('/file', methods=['GET'])
def handle_file():
    email = get_auth().auth()['email']
    database_api = get_database()
    dataset_id = request.args.get('id', '')
    url = request.args.get('file')
    if dataset_id == '':  # allow if file is in static directory
        static_dirs = os.environ.get(CIRRO_STATIC_DIR)
        if static_dirs is None:
            return 'Not authorized', 401
        static_dirs = static_dirs.split(',')
        if os.path.dirname(url) in static_dirs:
            return send_file(url)
        return 'Not authorized', 401
    else:
        dataset = database_api.get_dataset(email, dataset_id)
        file_path = get_file_path(url, dataset['url'])
        file_path = map_url(file_path)
    return send_file(file_path)


@cirro_blueprint.route('/user', methods=['GET'])
def handle_user():
    email = get_auth().auth()['email']
    database_api = get_database()
    user = database_api.user(email)
    return json_response(user)


@cirro_blueprint.route('/data', methods=['POST'])
def handle_data():
    json_request = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(json_request)
    return json_response(
        data_processing.handle_data(dataset_api=dataset_api,
                                    dataset=dataset,
                                    embedding_list=json_request.get('embedding'),
                                    values=json_request.get('values'),
                                    stats=json_request.get('stats'),
                                    grouped_stats=json_request.get('groupedStats'),
                                    selection=json_request.get('selection')))


@cirro_blueprint.route('/selection', methods=['POST'])
def handle_selection():
    # selection includes stats, coordinates
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    data_filter = content.get('filter')
    return json_response(
        data_processing.handle_selection(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter,
                                         measures=content.get('measures', []),
                                         dimensions=content.get('dimensions', []),
                                         embeddings=content.get('embeddings', [])))


@cirro_blueprint.route('/selected_ids', methods=['POST'])
def handle_selected_ids():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    data_filter = content.get('filter')
    return json_response(
        data_processing.handle_selection_ids(dataset_api=dataset_api, dataset=dataset, data_filter=data_filter))


# List available datasets
@cirro_blueprint.route('/datasets', methods=['GET'])
def handle_datasets():
    email = get_auth().auth()['email']
    database_api = get_database()
    return json_response(database_api.datasets(email))


@cirro_blueprint.route('/dataset', methods=['GET', 'POST', 'PUT', 'DELETE'])
def handle_dataset():
    email = get_auth().auth()['email']
    database_api = get_database()
    # POST=new dataset, PUT=update dataset, DELETE=delete, GET=get dataset info
    if request.method == 'PUT' or request.method == 'POST':
        if request.content_type == 'application/json':
            d = request.get_json(force=True, cache=False)
        else:
            d = dict()
            for key in request.form:
                d[key] = request.form[key]

        dataset_id = d.get('id')
        dataset_name = d.get('name')
        url = d.get('url')  # e.g. gs://foo/a/b/
        readers = d.pop('readers') if 'readers' in d else None
        file = None
        if request.content_type != 'application/json':
            file = request.files.get('file')  # file upload
        if readers is not None and not isinstance(readers, list):
            import json
            readers = json.loads(readers)
        if request.method == 'PUT' and dataset_id is None:  # update
            return 'Please supply an id', 400
        if request.method == 'POST' and dataset_name == '':  # new
            return 'Must supply dataset name', 400

        if url is not None and os.environ.get(CIRRO_UPLOAD) is not None:
            url = copy_url(url)
            d['url'] = url

        if file is not None:
            if os.environ.get(CIRRO_UPLOAD) is None:
                return 'Upload not supported', 400
            url = upload_file(file)
            d['url'] = url
        if request.method == 'POST' and url is None:  # new
            return 'Must supply dataset URL', 400
        dataset_id = database_api.upsert_dataset(email=email, readers=readers, dataset=d)
        return json_response({'id': dataset_id})
    elif request.method == 'DELETE':
        content = request.get_json(force=True, cache=False)
        dataset_id = content.get('id', '')
        database_api.delete_dataset(email, dataset_id)
        return json_response('', 204)
    elif request.method == 'GET':
        dataset_id = request.args.get('id', '')
        return json_response(database_api.get_dataset(email, dataset_id, True))


@cirro_blueprint.route('/module', methods=['POST'])
def handle_module_score():
    content = request.get_json(cache=False)
    email, dataset = get_email_and_dataset(content)
    features = content.get('features')
    adata = dataset_api.read_dataset(dataset, keys=dict(X=features))
    return json_response(adata.X.mean(axis=1))


@cirro_blueprint.route('/job', methods=['GET', 'DELETE', 'POST'])
def handle_job():
    email = get_auth().auth()['email']
    database_api = get_database()
    if request.method == 'DELETE':
        content = request.get_json(force=True, cache=False)
        job_id = content.get('id', '')
        database_api.delete_job(email, job_id)
        delete_job(job_id)
        return json_response('', 204)
    elif request.method == 'POST':
        if os.environ.get('GAE_APPLICATION') is None:  # TODO
            content = request.get_json(force=True, cache=False)
            email, dataset = get_email_and_dataset(content)
            params = content.get('params')
            job_type = content.get('type')
            job_name = content.get('name')
            return dict(id=submit_job(database_api=database_api, dataset_api=dataset_api, email=email, dataset=dataset,
                                      job_name=job_name, job_type=job_type, params=params))
        else:
            raise ValueError('Submit job not supported on GAE')
    else:
        job_id = request.args['id']
        c = request.args['c']
        is_precomputed = job_id.startswith('cirro-')
        if c == 'status' or c == 'params':
            if is_precomputed:
                job = dict(status='complete') if c == 'status' else dict()
            else:
                job = database_api.get_job(email=email, job_id=job_id, return_type=c)
            if job is None:
                return json_response('', 404)  # job deleted
            return json_response(job, 200)
        if c != 'result':
            raise ValueError('c must be one of status, params, or result')
        if is_precomputed:  # precomputed result
            email = get_auth().auth()['email']
            suggested_dataset_id = request.args['ds']
            dataset = database_api.get_dataset(email, suggested_dataset_id)
            # precomputed results need to be a child of dataset
            dataset['url'] = map_url(dataset['url'])
            job_result = dataset_api.get_result(dataset, job_id)
            if get_scheme(job_result) == 'file' and not os.path.exists(job_result):
                return Response(job_result, content_type='application/json')
            else:
                return send_file(job_result)
        job = database_api.get_job(email=email, job_id=job_id, return_type=c)
        if job is None:
            return json_response('', 404)  # job deleted
        import anndata
        if isinstance(job, dict) and 'url' in job:
            url = job['url']
            content_type = job.get('content-type')
            if content_type == 'application/h5ad' or content_type == 'application/zarr':
                if content_type == 'application/h5ad':
                    with get_fs(url).open(url, mode='rb') as f:
                        adata = anndata.read(f)
                else:
                    adata = anndata.read_zarr(get_fs(url).get_mapper(url))
                adata_df = adata_to_df(adata)
                return Response(adata_df.to_json(double_precision=2, orient='records'), content_type='application/json')
            else:
                # URL to JSON or text
                return send_file(url)
        elif isinstance(job, dict):
            return json_response(job)
        elif isinstance(job, anndata.AnnData):
            return Response(adata_to_df(job).to_json(double_precision=2, orient='records'),
                            content_type='application/json')
        return job


@cirro_blueprint.route('/jobs', methods=['GET'])
def handle_jobs():
    email = get_auth().auth()['email']
    database_api = get_database()
    ds_id = request.args.get('id', '')
    return json_response(database_api.get_jobs(email=email, dataset_id=ds_id))

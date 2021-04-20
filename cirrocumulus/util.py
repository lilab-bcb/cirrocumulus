# import pandas.io.json as json
import pandas._libs.json as ujson
from flask import make_response

from .file_system_adapter import get_scheme


def to_json(data, response=200):
    # response = make_response(simplejson.dumps(data, check_circular=True), response)
    # response = make_response(json.dumps(data), response)
    s = ujson.dumps(data, double_precision=2, orient='values')
    # s = nujson.dumps(data, double_precision=1)
    response = make_response(s, response)
    response.headers['Content-Type'] = 'application/json'
    return response


def get_email_domain(email):
    at_index = email.find('@')
    domain = None
    if at_index != -1:
        domain = email[at_index + 1:]
    return domain


def load_dataset_schema(url):
    import os

    import fsspec
    import json

    def get_extension(path):
        name, ext = os.path.splitext(path)
        if ext == '.gz':
            name, ext = os.path.splitext(name)
            if ext == '.json':
                ext = '.json.gz'
        return ext

    extension = get_extension(url)
    json_schema = None
    if extension in ['.json', '.json.gz', '']:
        scheme = get_scheme(url)
        fs = fsspec.filesystem(scheme)
        if extension == '':
            url = os.path.join(url, 'index.json.gz')
            extension = get_extension(url)
        if extension == '.json.gz':
            import gzip
            with gzip.open(fs.open(url)) as f:
                json_schema = json.load(f)
        else:
            with fs.open(url) as f:
                json_schema = json.load(f)
    return json_schema

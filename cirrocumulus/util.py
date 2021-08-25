# import pandas.io.json as json
from urllib.parse import urlparse

import fsspec
import pandas._libs.json as ujson
from flask import make_response


def get_scheme(path):
    pr = urlparse(path)
    if len(pr.scheme) <= 1:  # for file paths: /foo/bar/test.h5ad or C:/foo/bar/test.h5ad
        return 'file'
    return pr.scheme


scheme_to_fs = {}


def get_fs(path):
    scheme = get_scheme(path)
    fs = scheme_to_fs.get(scheme, None)
    if fs is not None:
        return fs
    fs = fsspec.filesystem(scheme)
    scheme_to_fs[scheme] = fs
    return fs


def to_json(data):
    return ujson.dumps(data, double_precision=2, orient='values')


def json_response(data, response=200):
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

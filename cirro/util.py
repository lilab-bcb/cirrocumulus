import gzip
import io

import ujson


def to_gzip(data):
    gzip_buffer = io.BytesIO()
    gzip_file = gzip.GzipFile(mode='wb', compresslevel=6, fileobj=gzip_buffer)
    gzip_file.write(data)
    gzip_file.close()
    return gzip_buffer.getvalue()


def to_json(data):
    return to_gzip(
        ujson.dumps(data, escape_forward_slashes=False, double_precision=1, ensure_ascii=False).encode('UTF-8'))


def get_json_headers():
    return {
            'Content-Type': 'application/json',
            'Content-Encoding': 'gzip',
            'Access-Control-Allow-Origin': '*'
    }


def handle_options(methods=['GET']):
    headers = {
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': ', '.join(methods),
            'Access-Control-Allow-Headers': 'Authorization',
            'Access-Control-Max-Age': '86400'
    }
    return '', 204, headers

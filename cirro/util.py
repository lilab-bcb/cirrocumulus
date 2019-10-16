import ujson
from flask import make_response


# def to_json(data, response=200):
#     response = make_response(simplejson.dumps(data, check_circular=False), response)
#     response.headers['Content-Type'] = 'application/json'
#     return response


def to_json(data, response=200):
    response = make_response(
        ujson.dumps(data, escape_forward_slashes=False, double_precision=1, ensure_ascii=False).encode('UTF-8'),
        response)
    response.headers['Content-Type'] = 'application/json'
    return response

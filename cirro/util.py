# import pandas.io.json as json
import pandas._libs.json as ujson
from flask import make_response


def to_json(data, response=200):
    # response = make_response(simplejson.dumps(data, check_circular=True), response)
    # response = make_response(json.dumps(data), response)
    s = ujson.dumps(data, double_precision=2, orient='values')
    # s = nujson.dumps(data, double_precision=1)
    response = make_response(s, response)
    response.headers['Content-Type'] = 'application/json'
    return response

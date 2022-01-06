import os

import cachecontrol
import google.auth.transport.requests
import requests
from flask import request
from google.oauth2 import id_token

from cirrocumulus.envir import CIRRO_AUTH_CLIENT_ID


class GoogleAuth:

    def __init__(self):
        session = requests.session()
        self.client_id = os.environ.get(CIRRO_AUTH_CLIENT_ID)
        cached_session = cachecontrol.CacheControl(session)
        self.cached_request = google.auth.transport.requests.Request(session=cached_session)

    def auth(self):
        token = request.headers.get('Authorization')
        if token is not None:
            if not token.startswith('Bearer '):
                raise ValueError('Token should start with Bearer ')
            token = token.split('Bearer ')[1]
        else:
            token = request.args.get('access_token')
        idinfo = id_token.verify_oauth2_token(token, self.cached_request, self.client_id)
        if idinfo['aud'] != self.client_id:
            raise ValueError('Wrong aud')
        if idinfo['iss'] not in ['accounts.google.com', 'https://accounts.google.com']:
            raise ValueError('Wrong issuer.')
        return idinfo

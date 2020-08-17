import cachecontrol
import google.auth.transport.requests
import requests
from flask import request
from google.oauth2 import id_token


class GoogleAuth:

    def __init__(self, clientId):
        session = requests.session()
        self.clientId = clientId
        cached_session = cachecontrol.CacheControl(session)
        self.cached_request = google.auth.transport.requests.Request(session=cached_session)

    @property
    def client_id(self):
        return self.clientId

    def auth(self):
        token = request.headers.get('Authorization')
        if token is not None:
            if not token.startswith('Bearer '):
                raise ValueError('Token should start with Bearer ')
            token = token.split('Bearer ')[1]
        else:
            token = request.args.get('access_token')
        idinfo = id_token.verify_oauth2_token(token, self.cached_request, self.clientId)
        if idinfo['aud'] != self.clientId:
            raise ValueError('Wrong aud')
        if idinfo['iss'] not in ['accounts.google.com', 'https://accounts.google.com']:
            raise ValueError('Wrong issuer.')
        return idinfo

import json

import cachecontrol
import google.auth.transport.requests
import requests
from flask import request
from google.oauth2 import id_token


class GoogleAuth:

    def __init__(self):
        with open('cirrocumulus/config.json', 'r') as f:
            self.config = json.load(f)
        session = requests.session()
        cached_session = cachecontrol.CacheControl(session)
        self.cached_request = google.auth.transport.requests.Request(session=cached_session)

    @property
    def client_id(self):
        return self.config['clientId']

    def auth(self):
        token = request.headers.get('Authorization')
        if not token.startswith('Bearer '):
            raise ValueError('Token should start with Bearer ')
        token = token.split('Bearer ')[1]
        idinfo = id_token.verify_oauth2_token(token, self.cached_request, self.config['clientId'])
        if idinfo['aud'] != self.config['clientId']:
            raise ValueError('Wrong aud')
        if idinfo['iss'] not in ['accounts.google.com', 'https://accounts.google.com']:
            raise ValueError('Wrong issuer.')
        return idinfo

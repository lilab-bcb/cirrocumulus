import os
import asyncio

from flask import request
from okta_jwt_verifier import BaseJWTVerifier, JWTUtils

from cirrocumulus.auth_exception import AuthException
from cirrocumulus.envir import CIRRO_AUTH_CLIENT_ID, CIRRO_AUTH_ISSUER


class OktaAuth:
    def __init__(self):
        self.loop = asyncio.get_event_loop()
        self.jwt_verifier = BaseJWTVerifier(
            issuer=os.environ.get(CIRRO_AUTH_ISSUER),
            client_id=os.environ.get(CIRRO_AUTH_CLIENT_ID),
            audience="api://default",
        )

    def auth(self):
        token = request.headers.get("Authorization")
        if token is not None:
            if not token.startswith("Bearer "):
                raise ValueError("Token should start with Bearer ")
            token = token.split("Bearer ")[1]
        else:
            token = request.args.get("access_token")
        try:
            result = JWTUtils.parse_token(token)
            self.loop.run_until_complete(self.jwt_verifier.verify_id_token(token))
            return result[1]
        except Exception:
            raise AuthException()

import werkzeug


class AuthException(werkzeug.exceptions.HTTPException):
    code = 401
    description = "Please authenticate"

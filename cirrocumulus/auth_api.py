class AuthAPI:

    def __init__(self):
        self.provider = None

    @property
    def client_id(self):
        return self.provider.client_id

    def auth(self):
        return self.provider.auth()

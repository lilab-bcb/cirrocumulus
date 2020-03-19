class NoAuth:

    @property
    def client_id(self):
        return ''

    def auth(self):
        return {'email': ''}

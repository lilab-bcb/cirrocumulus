class LocalDbAPI:

    def __init__(self, paths):
        self.paths = paths

    def server(self):
        return {}

    def user(self, email):
        return {}

    def datasets(self, email):
        results = []
        for path in self.paths:
            results.append({'id': path, 'name': path})
        return results

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        return {'id': dataset_id, 'name': dataset_id, 'url': dataset_id}

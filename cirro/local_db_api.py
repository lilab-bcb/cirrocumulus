import os

from cirro.entity import Entity


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
            results.append({'id': path, 'name': os.path.splitext(os.path.basename(path))[0]})
        return results

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        result = Entity(dataset_id, {'name': os.path.splitext(os.path.basename(dataset_id))[0], 'url': dataset_id})
        return result

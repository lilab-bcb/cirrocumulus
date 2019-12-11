import json
import os

from cirro.entity import Entity


class LocalDbAPI:

    def __init__(self, paths):
        self.paths = paths

    def server(self):
        return dict(canWrite=True)

    def user(self, email):
        return {}

    def create_dataset_meta(self, path):
        result = {'id': path, 'url': path, 'name': os.path.splitext(os.path.basename(path))[0]}
        json_file = os.path.join(path, 'index.json')
        if os.path.isdir(path) and os.path.exists(json_file):
            result['url'] = json_file
            with open(json_file, 'rt') as f:
                result.update(json.load(f))
        return result

    def dataset_filters(self, email, dataset_id):
        return []

    def datasets(self, email):
        results = []
        for path in self.paths:
            results.append(self.create_dataset_meta(path))
        return results

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        result = Entity(dataset_id, self.create_dataset_meta(dataset_id))
        return result

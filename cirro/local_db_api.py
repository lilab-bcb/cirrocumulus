import json
import os

from cirro.entity import Entity


class LocalDbAPI:

    def __init__(self, paths):
        self.paths = paths
        self.dataset_filter = {}
        self.counter = 0

    def server(self):
        return dict(canWrite=True)

    def user(self, email):
        return {}

    def create_dataset_meta(self, path):
        result = {'id': path, 'url': path, 'name': os.path.splitext(os.path.basename(path))[0]}
        if os.path.basename(path).endswith('.json'):
            with open(path, 'rt') as f:
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

    def dataset_filters(self, email, dataset_id):
        results = []
        for key in self.dataset_filter:
            results.append({'id': key, 'name': self.dataset_filter[key]['name']})
        return results

    def delete_dataset_filter(self, email, filter_id):
        del self.dataset_filter[filter_id]

    def get_dataset_filter(self, email, filter_id):
        return self.dataset_filter[filter_id]

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):

        if filter_id is None:
            self.counter += 1
            filter_id = str(self.counter)

        entity = self.dataset_filter.get(filter_id)
        if entity is None:
            entity = {}
            self.dataset_filter[filter_id] = entity
        if filter_name is not None:
            entity['name'] = filter_name
        if dataset_filter is not None:
            entity['value'] = json.dumps(dataset_filter)
        if email is not None:
            entity['email'] = email
        if dataset_id is not None:
            entity['dataset_id'] = dataset_id
        if filter_notes is not None:
            entity['notes'] = filter_notes
        return filter_id

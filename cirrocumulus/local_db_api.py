import json
import os

from cirrocumulus.entity import Entity


def create_dataset_meta(path):
    result = {'id': path, 'url': path, 'name': os.path.splitext(os.path.basename(path))[0]}
    if os.path.basename(path).endswith('.json'):
        with open(path, 'rt') as f:
            result.update(json.load(f))
    return result


class LocalDbAPI:

    def __init__(self, path):
        self.path = path
        self.dataset_filter_path = os.path.splitext(path)[0] + '_filters.json'
        self.dataset_filter = {}
        if os.path.exists(self.dataset_filter_path) and os.path.getsize(self.dataset_filter_path) > 0:
            with open(self.dataset_filter_path, 'rt') as f:
                self.dataset_filter.update(json.load(f))
        self.meta = create_dataset_meta(self.path)

    def server(self):
        return dict(canWrite=True)

    def user(self, email):
        return {}

    def datasets(self, email):
        results = []
        results.append(self.meta)
        return results

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        result = Entity(dataset_id, self.meta)
        return result

    def dataset_filters(self, email, dataset_id):
        results = []
        for key in self.dataset_filter:
            r = self.dataset_filter[key]
            r['id'] = key
            results.append(r)
        return results

    def delete_dataset_filter(self, email, filter_id):
        del self.dataset_filter[filter_id]
        self.__write_dataset_filter()

    def get_dataset_filter(self, email, filter_id):
        return self.dataset_filter[filter_id]

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        if filter_id is None:
            import uuid
            filter_id = str(uuid.uuid4())

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
        self.__write_dataset_filter()
        return filter_id

    def __write_dataset_filter(self):
        with open(self.dataset_filter_path, 'wt') as f:
            json.dump(self.dataset_filter, f)

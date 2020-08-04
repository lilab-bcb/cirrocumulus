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
        self.json_data = {}
        basename = os.path.splitext(path)[0]
        old_path = basename + '_filters.json'
        self.json_path = basename + '.json'
        if os.path.exists(old_path) and os.path.getsize(old_path) > 0:
            with open(old_path, 'rt') as f:
                self.json_data['filters'] = json.load(f)

        if os.path.exists(self.json_path) and os.path.getsize(self.json_path) > 0:
            with open(self.json_path, 'rt') as f:
                self.json_data.update(json.load(f))
        self.meta = create_dataset_meta(self.path)
        if 'filters' not in self.json_data:
            self.json_data['filters'] = {}
        if 'categories' not in self.json_data:
            self.json_data['categories'] = {}

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

    def category_names(self, dataset_id):
        results = []
        categories = self.json_data['categories']
        for category_name in categories:
            category = categories[category_name]
            for category_key in category:
                r = dict(category=category_name, original=category_key, new=category[category_key]['new'])
                results.append(r)
        return results

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        category_entity = self.json_data['categories'].get(category)
        if category_entity is None:
            category_entity = {}
            self.json_data['categories'][category] = category_entity
        if new_name == '':
            if original_name in category_entity:
                del category_entity[original_name]
        else:
            entity = dict(new=new_name)
            category_entity[original_name] = entity
            if dataset_id is not None:
                entity['dataset_id'] = dataset_id
            if email is not None:
                entity['email'] = email
        self.__write_json()

    def dataset_filters(self, email, dataset_id):
        results = []
        filters = self.json_data['filters']
        for key in filters:
            r = filters[key]
            r['id'] = key
            results.append(r)
        return results

    def delete_dataset_filter(self, email, filter_id):
        del self.json_data['filters'][filter_id]
        self.__write_json()

    def get_dataset_filter(self, email, filter_id):
        return self.json_data['filters'][filter_id]


    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        if filter_id is None:
            import uuid
            filter_id = str(uuid.uuid4())

        entity = self.json_data['filters'].get(filter_id)
        if entity is None:
            entity = {}
            self.json_data['filters'][filter_id] = entity
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
        self.__write_json()
        return filter_id

    def __write_json(self):
        with open(self.json_path, 'wt') as f:
            json.dump(self.json_data, f)

import datetime
import json
from urllib.parse import urlparse

import fsspec
from google.cloud import datastore

from .invalid_usage import InvalidUsage

DATASET = 'Dataset'
DATASET_FILTER = 'Dataset_Filter'
USER = 'User'


class FirestoreDatastore:

    def __init__(self):
        self.datastore_client = datastore.Client()

    def server(self):
        client = self.datastore_client
        server = dict(email=client.project + '@appspot.gserviceaccount.com')
        return server

    def user(self, email):
        client = self.datastore_client
        key = client.key(USER, email)
        user = client.get(key)
        if user is None:
            user = datastore.Entity(client.key(USER, email))
        user.update({'last_login': datetime.datetime.now()})
        client.put(user)
        return user

    def datasets(self, email):
        client = self.datastore_client
        query = client.query(kind=DATASET)
        query.add_filter('readers', '=', email)

        results = []
        for result in query.fetch():
            results.append({'id': result.id, 'name': result['name'], 'owner': email in result['owners']})
        return results

    def dataset_filters(self, email, dataset_id):
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        query = client.query(kind=DATASET_FILTER)
        query.add_filter('dataset_id', '=', int(dataset_id))
        results = []
        for result in query.fetch():
            results.append(result)
        return results

    def delete_dataset_filter(self, email, filter_id):
        client = self.datastore_client
        key = client.key(DATASET_FILTER, int(filter_id))
        dataset_filter = client.get(key)
        self.__get_key_and_dataset(email, dataset_filter['dataset_id'])
        client.delete(key)

    def get_dataset_filter(self, email, filter_id):
        client = self.datastore_client
        dataset_filter = client.get(client.key(DATASET_FILTER, int(filter_id)))
        self.__get_key_and_dataset(email, dataset_filter['dataset_id'])
        return dataset_filter

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        dataset_id = int(dataset_id)
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        if filter_id is None:
            dataset_filter_entity = datastore.Entity(client.key(DATASET_FILTER),
                exclude_from_indexes=['value'])
        else:
            filter_id = int(filter_id)
            key = client.key(DATASET_FILTER, filter_id)
            dataset_filter_entity = client.get(key)
        entity_update = {}
        if filter_name is not None:
            entity_update['name'] = filter_name
        if dataset_filter is not None:
            entity_update['value'] = json.dumps(dataset_filter)
        if email is not None:
            entity_update['email'] = email
        if dataset_id is not None:
            entity_update['dataset_id'] = dataset_id
        if filter_notes is not None:
            entity_update['notes'] = filter_notes

        dataset_filter_entity.update(entity_update)
        client.put(dataset_filter_entity)
        return dataset_filter_entity.id

    def __get_key_and_dataset(self, email, dataset_id, ensure_owner=False):
        client = self.datastore_client
        key = client.key(DATASET, int(dataset_id))
        dataset = client.get(key)
        if dataset is None:
            raise InvalidUsage('Please provide a valid id', 400)
        readers = dataset.get('readers')
        if email not in readers:
            raise InvalidUsage('Not authorized', 403)
        if ensure_owner and email not in dataset['owners']:
            raise InvalidUsage('Not authorized', 403)
        return key, dataset

    def delete_dataset(self, email, dataset_id):
        client = self.datastore_client
        key, dataset = self.__get_key_and_dataset(email, dataset_id, True)
        client.delete(key)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        key, dataset = self.__get_key_and_dataset(email, dataset_id, ensure_owner)
        return dataset

    def upsert_dataset(self, email, dataset_id, dataset_name, url, readers):
        client = self.datastore_client
        if dataset_id is not None:
            key, dataset = self.__get_key_and_dataset(email, dataset_id, True)
        else:
            dataset = datastore.Entity(client.key(DATASET), exclude_from_indexes=['url', 'precomputed'])
        readers = set(readers)
        if email in readers:
            readers.remove(email)
        readers.add(email)
        update_dict = {'name': dataset_name,
                       'readers': list(readers),
                       'url': url}
        if url.endswith('.json'):
            pr = urlparse(url)
            fs = fsspec.filesystem(pr.scheme if not pr.scheme == '' else 'file')
            with fs.open(url) as s:
                s = json.load(s)
                update_dict['precomputed'] = s.get('precomputed', False)
                update_dict['nObs'] = s.get('nObs')

        if dataset_id is None:  # new dataset
            update_dict['owners'] = [email]

        dataset.update(update_dict)
        client.put(dataset)
        dataset_id = dataset.id
        return dataset_id

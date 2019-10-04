import datetime

from google.cloud import datastore

from .invalid_usage import InvalidUsage

DATASET = 'Dataset'
USER = 'User'


class DatastoreAPI:

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
        query.add_filter('roles.{}'.format(email), '>', '')
        results = []
        for result in query.fetch():
            results.append({'id': result.id, 'name': result['name'], 'owner': result['roles'][email] == 'owner'})
        return results

    def __get_key_dataset(self, email, dataset_id, ensure_owner):
        client = self.datastore_client
        key = client.key(DATASET, int(dataset_id))
        dataset = client.get(key)
        if dataset is None:
            raise InvalidUsage('Please provide a valid id', 400)
        roles = dataset.get('roles')
        if email not in roles:
            raise InvalidUsage('Not authorized', 403)
        if ensure_owner and roles.get(email, '') != 'owner':
            raise InvalidUsage('Not authorized', 403)
        return key, dataset

    def delete_dataset(self, email, dataset_id):
        client = self.datastore_client
        key, dataset = self.__get_key_dataset(email, dataset_id, True)
        client.delete(key)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        key, dataset = self.__get_key_dataset(email, dataset_id, ensure_owner)
        return dataset

    def upsert_dataset(self, email, dataset_id, dataset_name, url, readers):
        client = self.datastore_client
        if dataset_id is not None:
            key, dataset = self.__get_key_dataset(email, dataset_id, True)
        else:
            dataset = datastore.Entity(client.key(DATASET))
            # if request_util.dataset_writer_collection.document(email) is None:
            #     return 'Not authorized to create datasets', 403
        if email in readers:
            readers.remove(email)
        roles = {}
        for user in readers:
            roles[user] = 'reader'
        roles[email] = 'owner'

        dataset.update({
                'name': dataset_name,
                'roles': roles,
                'url': url
        })
        client.put(dataset)
        dataset_id = dataset.id
        return dataset_id

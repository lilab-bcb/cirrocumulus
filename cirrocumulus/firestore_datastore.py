import datetime
import json

from cirrocumulus.util import get_email_domain
from google.cloud import datastore

from .invalid_usage import InvalidUsage

DATASET = 'Dataset'
CAT_NAME = 'Cat_Name'
DATASET_FILTER = 'Dataset_Filter'
FEATURE_SET = 'F_Set'
USER = 'User'
JOB = 'Job'
JOB_RESULT = 'Job_Result'


class FirestoreDatastore:

    def __init__(self):
        self.datastore_client = datastore.Client()

    def __get_key_and_dataset(self, email, dataset_id, ensure_owner=False):
        client = self.datastore_client
        key = client.key(DATASET, int(dataset_id))
        dataset = client.get(key)
        if dataset is None:
            raise InvalidUsage('Please provide a valid id', 400)
        readers = dataset.get('readers')

        domain = get_email_domain(email)
        if email not in readers and domain not in readers:
            raise InvalidUsage('Not authorized', 403)
        if ensure_owner and email not in dataset['owners']:
            raise InvalidUsage('Not authorized', 403)
        return key, dataset

    def server(self):
        client = self.datastore_client
        return dict(mode='server', email=client.project + '@appspot.gserviceaccount.com')

    def category_names(self, email, dataset_id):
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        client = self.datastore_client
        query = client.query(kind=CAT_NAME)
        query.add_filter('dataset_id', '=', dataset_id)
        results = []
        for result in query.fetch():
            results.append(result)
        return results

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        client = self.datastore_client
        key = client.key(CAT_NAME, str(dataset_id) + '-' + str(category) + '-' + str(original_name))

        if new_name == '':
            client.delete(key)
        else:
            entity = datastore.Entity(key=key)
            entity.update(dict(category=category, dataset_id=dataset_id, original=original_name, new=new_name))
            client.put(entity)
            return entity.id

    def user(self, email):
        client = self.datastore_client
        key = client.key(USER, email)
        user = client.get(key)
        if user is None:
            user = datastore.Entity(client.key(USER, email))
        user.update({'last_login': datetime.datetime.utcnow()})
        client.put(user)
        return user

    def datasets(self, email):
        client = self.datastore_client
        query = client.query(kind=DATASET)
        query.add_filter('readers', '=', email)
        results = []
        for result in query.fetch():
            results.append(
                {'id': result.id, 'name': result['name'], 'title': result.get('title'),
                 'owner': email in result['owners'], 'species': result.get('species'), 'url': result['url']})
        domain = get_email_domain(email)
        if domain is not None:
            query = client.query(kind=DATASET)
            query.add_filter('readers', '=', domain)
            unique_ids = set()
            for result in results:
                unique_ids.add(result['id'])
            for result in query.fetch():
                if result.id not in unique_ids:
                    results.append({'id': result.id, 'name': result['name'], 'title': result.get('title'),
                                    'species': result.get('species'),
                                    'owner': email in result['owners'],
                                    'url': result['url']})
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

    def delete_dataset_filter(self, email, dataset_id, filter_id):
        client = self.datastore_client
        key = client.key(DATASET_FILTER, int(filter_id))
        dataset_filter = client.get(key)
        self.__get_key_and_dataset(email, dataset_filter['dataset_id'])
        client.delete(key)

    def get_dataset_filter(self, email, dataset_id, filter_id):
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

    def delete_dataset(self, email, dataset_id):
        client = self.datastore_client
        key, dataset = self.__get_key_and_dataset(email, dataset_id, True)
        client.delete(key)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        _, dataset = self.__get_key_and_dataset(email, dataset_id, ensure_owner)
        dataset['id'] = dataset.id
        return dataset

    def is_importer(self, email):
        # TODO check if user can modify dataset
        client = self.datastore_client
        user = client.get(client.key(USER, email))
        if 'importer' not in user:
            raise False
        return user['importer']

    def upsert_dataset(self, email, dataset_id, dataset_name=None, url=None, readers=None, description=None, title=None,
                       species=None):
        client = self.datastore_client
        update_dict = {}
        if readers is not None:
            readers = set(readers)
            if email in readers:
                readers.remove(email)
            readers.add(email)
            update_dict['readers'] = list(readers)

        if dataset_id is not None:  # only owner can update
            key, dataset = self.__get_key_and_dataset(email, dataset_id, True)
        else:
            dataset = datastore.Entity(client.key(DATASET), exclude_from_indexes=['url'])
            user = client.get(client.key(USER, email))
            if 'importer' not in user or not user['importer']:
                raise InvalidUsage('Not authorized', 403)
            update_dict['owners'] = [email]
            if 'readers' not in update_dict:
                update_dict['readers'] = [email]
        if dataset_name is not None:
            update_dict['name'] = dataset_name
        if url is not None:
            update_dict['url'] = url

        if description is not None:
            update_dict['description'] = description

        if title is not None:
            update_dict['title'] = title
        if species is not None:
            update_dict['species'] = species
        dataset.update(update_dict)
        client.put(dataset)
        dataset_id = dataset.id
        return dataset_id

    def get_feature_sets(self, email, dataset_id):
        client = self.datastore_client
        query = client.query(kind=FEATURE_SET)
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        query.add_filter('dataset_id', '=', dataset_id)
        results = []
        for result in query.fetch():
            results.append(
                dict(id=result.id, category=result['category'], name=result['name'], features=result['features']))
        return results

    def delete_feature_set(self, email, dataset_id, set_id):
        client = self.datastore_client
        key = client.key(FEATURE_SET, int(set_id))
        result = client.get(key)
        self.__get_key_and_dataset(email, result['dataset_id'])
        client.delete(key)

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        dataset_id = int(dataset_id)
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        if set_id is None:
            entity = datastore.Entity(client.key(FEATURE_SET))
        else:
            set_id = int(set_id)
            key = client.key(FEATURE_SET, set_id)
            entity = client.get(key)
        entity_update = {}
        if name is not None:
            entity_update['name'] = name
        if category is not None:
            entity_update['category'] = category
        if email is not None:
            entity_update['email'] = email
        if dataset_id is not None:
            entity_update['dataset_id'] = dataset_id
        if features is not None:
            entity_update['features'] = features

        entity.update(entity_update)
        client.put(entity)
        return entity.id

    def create_job(self, email, dataset_id, job_name, job_type, params):
        dataset_id = int(dataset_id)
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        entity = datastore.Entity(client.key(JOB), exclude_from_indexes=["params"])
        import json
        entity.update(dict(dataset_id=dataset_id, email=email, name=job_name, type=job_type, status='',
            submitted=datetime.datetime.utcnow(), params=json.dumps(params)))
        client.put(entity)
        job_id = entity.id

        entity = datastore.Entity(client.key(JOB_RESULT, job_id), exclude_from_indexes=["result"])
        entity.update(dict(dataset_id=dataset_id))
        client.put(entity)
        return str(job_id)

    def get_job(self, email, job_id, return_result):
        client = self.datastore_client
        job_id = int(job_id)
        key = client.key(JOB_RESULT if return_result else JOB, job_id)
        result = client.get(key)
        self.__get_key_and_dataset(email, result['dataset_id'])
        return dict(status=result['result' if return_result else 'status'])

    def update_job(self, email, job_id, status, result):
        client = self.datastore_client
        job_id = int(job_id)
        is_complete = result is not None
        key = client.key(JOB, job_id)
        entity = client.get(key)
        if not is_complete:
            self.__get_key_and_dataset(email, entity['dataset_id'])
        entity.update(dict(status=status))
        client.put(entity)

        if is_complete:
            key = client.key(JOB_RESULT, job_id)
            entity = client.get(key)
            entity.update(dict(result=result))
            client.put(entity)

    def get_jobs(self, email, dataset_id):
        client = self.datastore_client
        query = client.query(kind=JOB)
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        query.add_filter('dataset_id', '=', dataset_id)
        results = []
        for result in query.fetch():
            results.append(
                dict(id=result.id, name=result['name'], type=result['type'], submitted=result['submitted'],
                    status=result['status'], email=result['email']))
        return results

    def delete_job(self, email, job_id):
        client = self.datastore_client
        job_id = int(job_id)
        key = client.key(JOB, job_id)
        result = client.get(key)
        if result['email'] == email:
            client.delete(key)
            client.delete(client.key(JOB_RESULT, job_id))
        else:
            raise InvalidUsage('Not authorized', 403)

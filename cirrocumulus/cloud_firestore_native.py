import datetime
import json
import os

from cirrocumulus.abstract_db import AbstractDB
from cirrocumulus.envir import SERVER_CAPABILITY_JOBS, CIRRO_EMAIL, CIRRO_JOB_RESULTS
from cirrocumulus.invalid_usage import InvalidUsage
from cirrocumulus.job_api import save_job_result_to_file
from cirrocumulus.util import get_email_domain
from google.cloud import datastore

DATASET = 'Dataset'
CAT_NAME = 'Cat_Name'
DATASET_FILTER = 'Dataset_Filter'
DATASET_VIEW = 'D_View'
FEATURE_SET = 'F_Set'
USER = 'User'
JOB = 'Job'


def get_datasets(results, email, query, unique_ids):
    for result in query.fetch():
        if result.id not in unique_ids:
            unique_ids.add(result.id)
            results.append(
                {'id': result.id, 'name': result['name'], 'title': result.get('title'),
                 'owner': email in result['owners'], 'species': result.get('species'), 'url': result['url']})


class CloudFireStoreNative(AbstractDB):

    def __init__(self):
        super().__init__()
        self.datastore_client = datastore.Client()
        os.environ[CIRRO_EMAIL] = self.datastore_client.project + '@appspot.gserviceaccount.com'

    def capabilities(self):
        c = super().capabilities()
        c[SERVER_CAPABILITY_JOBS] = False
        return c

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

    def __get_entity_list(self, email, dataset_id, kind, keys):
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        query = client.query(kind=kind)
        query.add_filter('dataset_id', '=', int(dataset_id))
        results = []
        for result in query.fetch():
            d = dict(id=result.id)
            for key in keys:
                d[key] = result.get(key)
            results.append(result)
        return results

    def __delete_entity(self, email, kind, entity_id):
        client = self.datastore_client
        key = client.key(kind, int(entity_id))
        e = client.get(key)
        self.__get_key_and_dataset(email, e['dataset_id'])
        client.delete(key)

    def __get_entity(self, email, kind, entity_id):
        client = self.datastore_client
        e = client.get(client.key(kind, int(entity_id)))
        self.__get_key_and_dataset(email, e['dataset_id'])
        return e

    def __upsert_entity(self, email, dataset_id, entity_id, kind, entity_update):
        dataset_id = int(dataset_id)
        client = self.datastore_client
        self.__get_key_and_dataset(email, dataset_id)
        if entity_id is None:
            e = datastore.Entity(client.key(kind),
                                 exclude_from_indexes=['value'])
        else:
            entity_id = int(entity_id)
            key = client.key(kind, entity_id)
            e = client.get(key)

        if email is not None:
            entity_update['email'] = email
        if dataset_id is not None:
            entity_update['dataset_id'] = dataset_id

        e.update(entity_update)
        client.put(e)
        return e.id

    def category_names(self, email, dataset_id):
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        client = self.datastore_client
        query = client.query(kind=CAT_NAME)
        query.add_filter('dataset_id', '=', dataset_id)
        results = {}
        for result in query.fetch():
            category = result['category']
            cat_results = results.get(category)
            if cat_results is None:
                cat_results = {}
                results[category] = cat_results
            cat_results[result['original']] = {'newValue': result.get('new'),
                                               'negativeMarkers': result.get('negativeMarkers'),
                                               'positiveMarkers': result.get('positiveMarkers'),
                                               'color': result.get('color')}

        return results

    def upsert_category_name(self, email, dataset_id, category, original_value, update):
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        client = self.datastore_client
        key = client.key(CAT_NAME, str(dataset_id) + '-' + str(category) + '-' + str(original_value))
        entity = datastore.Entity(key=key)
        update['category'] = category
        update['dataset_id'] = dataset_id
        update['original'] = original_value
        if 'newValue' in update:  # backwards compatibility
            update['new'] = update.pop('newValue')
        entity.update(update)
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
        unique_ids = set()
        get_datasets(results, email, query, unique_ids)
        domain = get_email_domain(email)
        if domain is not None:
            query = client.query(kind=DATASET)
            query.add_filter('readers', '=', domain)
            # TODO make this more efficient
            get_datasets(results, email, query, unique_ids)

            query = client.query(kind=DATASET)
            query.add_filter('readers', '=', 'allUsers')
            get_datasets(results, email, query, unique_ids)
        return results

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

    def upsert_dataset(self, email, readers, dataset):
        client = self.datastore_client

        if readers is not None:
            readers = set(readers)
            if email in readers:
                readers.remove(email)
            readers.add(email)
            dataset['readers'] = list(readers)
        dataset_id = dataset.get('id')
        if dataset_id is not None:  # only owner can update
            key, dataset_entity = self.__get_key_and_dataset(email, dataset_id, True)
        else:  # new dataset
            dataset_entity = datastore.Entity(client.key(DATASET), exclude_from_indexes=['url'])
            user = client.get(client.key(USER, email))
            if 'importer' not in user or not user['importer']:
                raise InvalidUsage('Not authorized', 403)
            dataset['owners'] = [email]

        dataset_entity.update(dataset)
        client.put(dataset_entity)
        dataset_id = dataset.id
        return dataset_id

    # feature sets
    def get_feature_sets(self, email, dataset_id):
        return self.__get_entity_list(email=email, dataset_id=dataset_id, kind=FEATURE_SET,
                                      keys=['category', 'name', 'features'])

    def delete_feature_set(self, email, dataset_id, set_id):
        return self.__delete_entity(email=email, kind=FEATURE_SET, entity_id=set_id)

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        entity_update = {}
        if name is not None:
            entity_update['name'] = name
        if category is not None:
            entity_update['category'] = category
        if features is not None:
            entity_update['features'] = features
        return self.__upsert_entity(email=email, dataset_id=dataset_id, entity_id=set_id, kind=FEATURE_SET,
                                    entity_update=entity_update)

    # views
    def dataset_views(self, email, dataset_id):
        return self.__get_entity_list(email=email, dataset_id=dataset_id, kind=DATASET_VIEW,
                                      keys=['name', 'notes', 'email', 'last_updated'])

    def delete_dataset_view(self, email, view_id):
        return self.__delete_entity(email=email, kind=DATASET_VIEW, entity_id=view_id)

    def get_dataset_view(self, email, view_id):
        return self.__get_entity(email=email, entity_id=view_id, kind=DATASET_VIEW)

    def upsert_dataset_view(self, email, dataset_id, view):
        view['last_updated'] = datetime.datetime.utcnow()
        if 'value' in view:
            view['value'] = json.dumps(view['value'])
        if email is not None:
            view['email'] = email
        if dataset_id is not None:
            view['dataset_id'] = dataset_id
        view_id = view.pop('id') if 'id' in view else None
        view_id = self.__upsert_entity(email=email, dataset_id=dataset_id, entity_id=view_id, kind=DATASET_VIEW,
                                       entity_update=view)
        return dict(id=view_id, last_updated=view['last_updated'])

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
        return str(job_id)

    def get_job(self, email, job_id, return_type):
        client = self.datastore_client
        job_id = int(job_id)
        key = client.key(JOB, job_id)
        result = client.get(key)
        self.__get_key_and_dataset(email, result['dataset_id'])
        if return_type == 'result':
            return result['result']
        elif return_type == 'status':
            return dict(status=result['status'])
        elif return_type == 'params':
            return dict(params=result.get('params'))

    def update_job(self, email, job_id, status, result):
        if not self.capabilities()[SERVER_CAPABILITY_JOBS]:
            return
        client = self.datastore_client
        job_id = int(job_id)
        is_complete = result is not None
        key = client.key(JOB, job_id)
        entity = client.get(key)
        if not is_complete:
            self.__get_key_and_dataset(email, entity['dataset_id'])

        if is_complete:
            if os.environ.get(CIRRO_JOB_RESULTS) is not None:  # save to directory
                result = save_job_result_to_file(result, job_id)
            else:
                raise ValueError('Please specify job result location')
            entity.update(dict(result=result, status=status))
            client.put(entity)
        else:
            entity.update(dict(status=status))
            client.put(entity)

    def get_jobs(self, email, dataset_id):
        client = self.datastore_client
        results = []

        query = client.query(kind=JOB)
        dataset_id = int(dataset_id)
        self.__get_key_and_dataset(email, dataset_id, False)
        query.add_filter('dataset_id', '=', dataset_id)

        for result in query.fetch():
            results.append(
                dict(id=result.id, name=result['name'], type=result['type'], submitted=result.get('submitted'),
                     status=result.get('status'), email=result.get('email')))
        return results

    def delete_job(self, email, job_id):
        client = self.datastore_client
        job_id = int(job_id)
        key = client.key(JOB, job_id)
        result = client.get(key)
        if result['email'] == email:
            client.delete(key)
        else:
            raise InvalidUsage('Not authorized', 403)

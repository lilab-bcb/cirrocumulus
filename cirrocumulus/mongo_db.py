import datetime
import os

import pandas._libs.json as ujson
from bson import ObjectId
from pymongo import MongoClient

from cirrocumulus.abstract_db import AbstractDB
from cirrocumulus.util import get_email_domain, get_fs
from .envir import CIRRO_DB_URI, CIRRO_AUTH_CLIENT_ID, CIRRO_JOB_RESULTS, SERVER_CAPABILITY_EDIT_DATASET, \
    SERVER_CAPABILITY_ADD_DATASET, SERVER_CAPABILITY_DELETE_DATASET, SERVER_CAPABILITY_LINKS, SERVER_CAPABILITY_JOBS, \
    SERVER_CAPABILITY_FEATURE_SETS, SERVER_CAPABILITY_RENAME_CATEGORIES
from .invalid_usage import InvalidUsage
from .job_api import save_job_result_to_file


def format_doc(d):
    d['id'] = str(d.pop('_id'))
    return d


class MongoDb(AbstractDB):

    def __init__(self):
        super().__init__()
        self.client = MongoClient(os.environ[CIRRO_DB_URI])
        self.db = self.client.get_default_database()
        self.db.categories.create_index('dataset_id')
        self.db.views.create_index('dataset_id')
        self.db.feature_sets.create_index('dataset_id')
        self.db.jobs.create_index('dataset_id')
        self.db.datasets.create_index('readers')
        self.fs = None

    def category_names(self, email, dataset_id):

        self.get_dataset(email, dataset_id)
        collection = self.db.categories
        results = {}
        for doc in collection.find(dict(dataset_id=dataset_id)):
            category = doc['category']
            cat_results = results.get(category)
            if cat_results is None:
                cat_results = {}
                results[category] = cat_results
            cat_results[doc['original']] = {'newValue': doc.get('new'),
                                            'negativeMarkers': doc.get('negativeMarkers'),
                                            'positiveMarkers': doc.get('positiveMarkers'),
                                            'color': doc.get('color')}
        return results

    def upsert_category_name(self, email, dataset_id, category, original_value, update):
        if not self.capabilities()[SERVER_CAPABILITY_RENAME_CATEGORIES]:
            return
        self.get_dataset(email, dataset_id)
        collection = self.db.categories
        key = str(dataset_id) + '-' + str(category) + '-' + str(original_value)
        update['category'] = category
        update['dataset_id'] = dataset_id
        update['original'] = original_value
        update['last_updated'] = datetime.datetime.utcnow()
        if 'newValue' in update:  # backwards compatibility
            update['new'] = update.pop('newValue')
        collection.update_one(dict(cat_id=key), {'$set': update}, upsert=True)

    def user(self, email):
        collection = self.db.users
        doc = collection.find_one(dict(email=email))
        if doc is None:
            collection.insert_one(dict(email=email))
            return {'id': email}
        else:
            return {'id': email,
                    'importer': doc.get('importer', False)}

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        collection = self.db.datasets
        auth_client_id = os.environ.get(CIRRO_AUTH_CLIENT_ID)

        if auth_client_id is None:  # allow unregistered URL
            try:
                dataset_id.index('://')
                return {
                    'id': dataset_id,
                    'name': dataset_id,
                    'url': dataset_id
                }
            except ValueError:
                pass
        doc = collection.find_one(dict(_id=ObjectId(dataset_id)))
        if doc is None:
            raise InvalidUsage('Please provide a valid id', 400)
        readers = doc.get('readers')
        domain = get_email_domain(email)
        if email not in readers and domain not in readers:
            raise InvalidUsage('Not authorized', 403)
        if ensure_owner and email not in doc['owners']:
            raise InvalidUsage('Not authorized', 403)
        doc['owner'] = 'owners' in doc and email in doc.pop('owners')
        doc['id'] = str(doc.pop('_id'))
        return doc

    def datasets(self, email):
        collection = self.db.datasets
        results = []
        domain = get_email_domain(email)
        if domain is None:
            query = dict(readers=email)
        else:
            query = dict(readers={'$in': [email, domain]})

        for doc in collection.find(query):
            doc['owner'] = 'owners' in doc and email in doc.pop('owners')
            doc['id'] = str(doc.pop('_id'))
            results.append(doc)
        return results

    # views
    def dataset_views(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.views
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id), dict(value=0)):
            results.append(format_doc(doc))
        return results

    def delete_dataset_view(self, email, view_id):
        if not self.capabilities()[SERVER_CAPABILITY_LINKS]:
            return
        collection = self.db.views
        doc = collection.find_one(dict(_id=ObjectId(view_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.delete_one(dict(_id=ObjectId(view_id)))

    def get_dataset_view(self, email, view_id):
        collection = self.db.views
        doc = collection.find_one(dict(_id=ObjectId(view_id)))
        self.get_dataset(email, doc['dataset_id'])
        return format_doc(doc)

    def upsert_dataset_view(self, email, dataset_id, view):
        if not self.capabilities()[SERVER_CAPABILITY_LINKS]:
            return
        self.get_dataset(email, dataset_id)
        collection = self.db.views
        view['last_updated'] = datetime.datetime.utcnow()

        if email is not None:
            view['email'] = email
        if dataset_id is not None:
            view['dataset_id'] = dataset_id
        if view.get('id') is None:
            return dict(id=str(collection.insert_one(view).inserted_id), last_updated=view['last_updated'])
        else:
            view_id = view.pop('id')
            collection.update_one(dict(_id=ObjectId(view_id)), {'$set': view})
            return dict(id=view_id, last_updated=view['last_updated'])

    def is_importer(self, email):
        # TODO check if user can modify dataset
        user = self.db.users.find_one(dict(email=email))
        if 'importer' not in user:
            raise False
        return user['importer']

    def delete_dataset(self, email, dataset_id):
        if not self.capabilities()[SERVER_CAPABILITY_DELETE_DATASET]:
            return

        self.get_dataset(email, dataset_id, True)
        self.db.datasets.delete_one(dict(_id=ObjectId(dataset_id)))
        self.db.filters.delete_many(dict(dataset_id=dataset_id))
        self.db.views.delete_many(dict(dataset_id=dataset_id))
        self.db.categories.delete_many(dict(dataset_id=dataset_id))
        self.db.feature_sets.delete_many(dict(dataset_id=dataset_id))

        for doc in self.db.jobs.find(dict(dataset_id=dataset_id), dict(result=1)):
            result = doc.get('result')
            self.delete_job_result(result)
        self.db.jobs.delete_many(dict(dataset_id=dataset_id))

    def upsert_dataset(self, email, readers, dataset):

        if dataset.get('id') is None and not self.capabilities()[SERVER_CAPABILITY_ADD_DATASET]:
            return
        if dataset.get('id') is not None and not self.capabilities()[SERVER_CAPABILITY_EDIT_DATASET]:
            return
        collection = self.db.datasets

        if readers is not None:
            readers = set(readers)
            if email in readers:
                readers.remove(email)
            readers.add(email)
            dataset['readers'] = list(readers)
        dataset['last_updated'] = datetime.datetime.utcnow()
        if dataset.get('id') is None:  # new dataset
            if email != '':
                user = self.db.users.find_one(dict(email=email))
                if 'importer' not in user or not user['importer']:
                    raise InvalidUsage('Not authorized', 403)
            dataset['owners'] = [email]
            if 'readers' not in dataset:
                dataset['readers'] = [email]
            return str(collection.insert_one(dataset).inserted_id)
        else:  # update
            self.get_dataset(email, dataset['id'], True)
            dataset_id = dataset.pop('id')
            collection.update_one(dict(_id=ObjectId(dataset_id)), {'$set': dataset})
            return dataset_id

    def get_feature_sets(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.feature_sets
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append(format_doc(doc))
        return results

    def delete_feature_set(self, email, dataset_id, set_id):
        if not self.capabilities()[SERVER_CAPABILITY_FEATURE_SETS]:
            return
        collection = self.db.feature_sets
        doc = collection.find_one(dict(_id=ObjectId(set_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.delete_one(dict(_id=ObjectId(set_id)))

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        if not self.capabilities()[SERVER_CAPABILITY_FEATURE_SETS]:
            return
        self.get_dataset(email, dataset_id)
        collection = self.db.feature_sets
        entity_update = {}
        entity_update['last_updated'] = datetime.datetime.utcnow()
        if name is not None:
            entity_update['name'] = name
        if features is not None:
            entity_update['features'] = features
        if email is not None:
            entity_update['email'] = email
        if dataset_id is not None:
            entity_update['dataset_id'] = dataset_id
        if category is not None:
            entity_update['category'] = category
        if set_id is None:
            return str(collection.insert_one(entity_update).inserted_id)
        else:
            collection.update_one(dict(_id=ObjectId(set_id)), {'$set': entity_update})
            return set_id

    def get_gridfs(self):
        import gridfs
        if self.fs is None:
            self.fs = gridfs.GridFS(self.db)
        return self.fs

    def create_job(self, email, dataset_id, job_name, job_type, params):
        if not self.capabilities()[SERVER_CAPABILITY_JOBS]:
            return
        self.get_dataset(email, dataset_id)
        collection = self.db.jobs
        job_id = str(collection.insert_one(
            dict(dataset_id=dataset_id, status='pending', name=job_name, email=email, type=job_type, params=params,
                 submitted=datetime.datetime.utcnow())).inserted_id)
        return job_id

    def get_job(self, email, job_id, return_type):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)))
        if doc is None:
            return None
        self.get_dataset(email, doc['dataset_id'])
        if return_type == 'result':
            result = doc.get('result')
            if isinstance(result, dict) and result.get('url') is not None:
                return result
            else:
                fs = self.get_gridfs()
                return fs.get(ObjectId(doc['result'])).read()
        elif return_type == 'status':
            return dict(status=doc['status'])
        elif return_type == 'params':
            return dict(params=doc.get('params'))

    def get_jobs(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.jobs
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id),
                                   dict(name=1, status=1, email=1, type=1, submitted=1, readonly=1)):
            results.append(
                dict(id=str(doc['_id']), name=doc['name'], status=doc['status'], type=doc['type'], email=doc['email'],
                     submitted=doc.get('submitted'), readonly=doc.get('readonly', False)))
        return results

    def annotate_job(self, email, job_id, annotations):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.update_one(dict(_id=ObjectId(job_id)),
                              {'$set': dict(annotations=annotations, last_updated=datetime.datetime.utcnow())})

    def update_job(self, email, job_id, status, result):
        if not self.capabilities()[SERVER_CAPABILITY_JOBS]:
            return
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)))
        self.get_dataset(email, doc['dataset_id'])
        if doc.get('readonly', False):
            raise InvalidUsage('Not authorized', 403)
        if result is not None:
            if os.environ.get(CIRRO_JOB_RESULTS) is not None:  # save to directory
                result = save_job_result_to_file(result, job_id)
            else:
                result = ujson.dumps(result, double_precision=2, orient='values')
                result = str(self.get_gridfs().put(result, encoding='ascii'))

        collection.update_one(dict(_id=ObjectId(job_id)), {'$set': dict(status=status, result=result)})

    def delete_job_result(self, result):
        if result is not None:
            if isinstance(result, dict) and result.get('url') is not None:
                fs = get_fs(result['url'])
                if fs.exists(result['url']):
                    fs.rm(result['url'], recursive=True)
            else:
                self.get_gridfs().delete(ObjectId(result))

    def delete_job(self, email, job_id):
        if not self.capabilities()[SERVER_CAPABILITY_JOBS]:
            return
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)), dict(email=1))
        if doc['email'] == email and not doc.get('readonly', False):
            collection.delete_one(dict(_id=ObjectId(job_id)))
            result = doc.get('result')
            self.delete_job_result(result)
        else:
            raise InvalidUsage('Not authorized', 403)

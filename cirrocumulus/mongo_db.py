import datetime
import json
import os

import pandas._libs.json as ujson
from bson import ObjectId
from cirrocumulus.abstract_db import AbstractDB
from cirrocumulus.util import get_email_domain, get_fs
from pymongo import MongoClient

from .envir import CIRRO_DB_URI, CIRRO_AUTH_CLIENT_ID, CIRRO_JOB_RESULTS
from .invalid_usage import InvalidUsage


def format_doc(d):
    d['id'] = str(d.pop('_id'))
    return d


class MongoDb(AbstractDB):

    def __init__(self):
        super().__init__()
        self.client = MongoClient(os.environ[CIRRO_DB_URI])
        self.db = self.client.get_default_database()

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
        self.get_dataset(email, dataset_id)
        collection = self.db.categories
        key = str(dataset_id) + '-' + str(category) + '-' + str(original_value)
        update['category'] = category
        update['dataset_id'] = dataset_id
        update['original'] = original_value
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

        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append(format_doc(doc))
        return results

    def delete_dataset_view(self, email, view_id):
        collection = self.db.views
        doc = collection.find_one(dict(_id=ObjectId(view_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.delete_one(dict(_id=ObjectId(view_id)))

    def get_dataset_view(self, email, view_id):
        collection = self.db.views
        doc = collection.find_one(dict(_id=ObjectId(view_id)))
        self.get_dataset(email, doc['dataset_id'])
        return format_doc(doc)

    def upsert_dataset_view(self, email, dataset_id, view_id, name, value):
        self.get_dataset(email, dataset_id)
        collection = self.db.views
        entity_update = {'last_updated': datetime.datetime.utcnow()}
        if name is not None:
            entity_update['name'] = name
        if value is not None:
            entity_update['value'] = json.dumps(value)
        if email is not None:
            entity_update['email'] = email
        if dataset_id is not None:
            entity_update['dataset_id'] = dataset_id
        if view_id is None:
            return str(collection.insert_one(entity_update).inserted_id)
        else:
            collection.update_one(dict(_id=ObjectId(view_id)), {'$set': entity_update})
            return view_id

    def delete_dataset(self, email, dataset_id):
        self.get_dataset(email, dataset_id, True)
        collection = self.db.datasets
        collection.delete_one(dict(_id=ObjectId(dataset_id)))
        self.db.filters.delete_many(dict(dataset_id=dataset_id))
        self.db.categories.delete_many(dict(dataset_id=dataset_id))

    def is_importer(self, email):
        # TODO check if user can modify dataset
        user = self.db.users.find_one(dict(email=email))
        if 'importer' not in user:
            raise False
        return user['importer']

    def upsert_dataset(self, email, dataset_id, dataset_name=None, url=None, readers=None, description=None, title=None,
                       species=None):
        collection = self.db.datasets
        update_dict = {}
        if dataset_name is not None:
            update_dict['name'] = dataset_name
        if url is not None:
            update_dict['url'] = url

        if readers is not None:
            readers = set(readers)
            if email in readers:
                readers.remove(email)
            readers.add(email)
            update_dict['readers'] = list(readers)
        if description is not None:
            update_dict['description'] = description
        if title is not None:
            update_dict['title'] = title
        if species is not None:
            update_dict['species'] = species

        if dataset_id is None:  # new dataset
            if email != '':
                user = self.db.users.find_one(dict(email=email))
                if 'importer' not in user or not user['importer']:
                    raise InvalidUsage('Not authorized', 403)
            update_dict['owners'] = [email]
            if 'readers' not in update_dict:
                update_dict['readers'] = [email]
            return str(collection.insert_one(update_dict).inserted_id)
        else:
            self.get_dataset(email, dataset_id, True)
            collection.update_one(dict(_id=ObjectId(dataset_id)), {'$set': update_dict})
            return dataset_id

    def get_feature_sets(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.feature_sets
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append(format_doc(doc))
        return results

    def delete_feature_set(self, email, dataset_id, set_id):
        collection = self.db.feature_sets
        doc = collection.find_one(dict(_id=ObjectId(set_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.delete_one(dict(_id=ObjectId(set_id)))

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        self.get_dataset(email, dataset_id)
        collection = self.db.feature_sets
        entity_update = {}
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

    def create_job(self, email, dataset_id, job_name, job_type, params):
        self.get_dataset(email, dataset_id)

        collection = self.db.jobs
        return str(collection.insert_one(
            dict(dataset_id=dataset_id, name=job_name, email=email, type=job_type, params=params,
                 submitted=datetime.datetime.utcnow())).inserted_id)

    def get_job(self, email, job_id, return_result):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)),
                                  {"result": 0} if not return_result else {'result': 1, "dataset_id": 1})
        self.get_dataset(email, doc['dataset_id'])
        if return_result:
            return doc['result']
        else:
            return dict(status=doc['status'])

    def get_jobs(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.jobs
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id), dict(name=1, status=1, email=1, type=1, submitted=1)):
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
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)))
        self.get_dataset(email, doc['dataset_id'])
        if doc.get('readonly', False):
            raise InvalidUsage('Not authorized', 403)
        if result is not None:

            if os.environ.get(CIRRO_JOB_RESULTS) is not None:  # save to directory
                url = os.path.join(os.environ[CIRRO_JOB_RESULTS], str(job_id) + '.json.gz')
                new_result = dict(url=url)
                new_result['content-type'] = result.pop('content-type')
                new_result['content-encoding'] = 'gzip'

                with get_fs(url).open(url, 'wt', compression='gzip') as out:
                    out.write(ujson.dumps(result, double_precision=2, orient='values'))
                result = new_result
            else:
                result = ujson.dumps(result, double_precision=2, orient='values')
        collection.update_one(dict(_id=ObjectId(job_id)), {'$set': dict(status=status, result=result)})

    def delete_job(self, email, job_id):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)), dict(email=1))
        if doc['email'] == email and not doc.get('readonly', False):
            collection.delete_one(dict(_id=ObjectId(job_id)))
            url = doc.get('url')
            if url is not None:
                get_fs(url).rm(url)
        else:
            raise InvalidUsage('Not authorized', 403)

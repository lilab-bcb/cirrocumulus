import json

from bson import ObjectId
from cirrocumulus.util import get_email_domain
from pymongo import MongoClient

from .invalid_usage import InvalidUsage


class MongoDb:

    def __init__(self, db_uri='mongodb://localhost:27017/', database='cirrocumulus', email=None):
        self.client = MongoClient(db_uri)
        self.db = self.client[database]
        self.email = email

    def server(self):
        d = dict(mode='server')
        if self.email is not None:
            d['email'] = self.email
        return d

    def category_names(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.categories
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append({'category': doc['category'], 'dataset_id': doc['dataset_id'], 'original': doc['original'],
                            'new': doc['new']})
        return results

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        self.get_dataset(email, dataset_id)
        collection = self.db.categories
        key = str(dataset_id) + '-' + str(category) + '-' + str(original_name)

        if new_name == '':
            collection.delete_one(dict(cat_id=key))
        else:
            collection.update_one(dict(cat_id=key),
                {'$set': dict(category=category, dataset_id=dataset_id, original=original_name, new=new_name)},
                upsert=True)
            return key

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
        doc = collection.find_one(dict(_id=ObjectId(dataset_id)))
        if doc is None:
            raise InvalidUsage('Please provide a valid id', 400)
        readers = doc.get('readers')
        domain = get_email_domain(email)
        if email not in readers and domain not in readers:
            raise InvalidUsage('Not authorized', 403)
        if ensure_owner and email not in doc['owners']:
            raise InvalidUsage('Not authorized', 403)
        return {
                'id': str(doc['_id']),
                'name': doc['name'],
                'readers': doc.get('readers'),
                'species': doc.get('species'),
                'description': doc.get('description'),
                'title': doc.get('title'),
                'url': doc['url'],
                'owner': 'owners' in doc and email in doc['owners']}

    def datasets(self, email):
        collection = self.db.datasets
        results = []
        domain = get_email_domain(email)
        if domain is None:
            query = dict(readers=email)
        else:
            query = dict(readers={'$in': [email, domain]})
        for doc in collection.find(query):
            results.append({'id': str(doc['_id']), 'name': doc['name'], 'title': doc.get('title'),
                            'owner': 'owners' in doc and email in doc['owners'], 'url': doc['url'],
                            'species': doc.get('species')})
        return results

    def dataset_filters(self, email, dataset_id):
        self.get_dataset(email, dataset_id)
        collection = self.db.filters
        results = []

        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append(
                {'id': str(doc['_id']), 'dataset_id': doc['dataset_id'], 'name': doc['name'], 'value': doc['value'],
                 'notes': doc['notes'], 'email': doc['email']})
        return results

    def delete_dataset_filter(self, email, dataset_id, filter_id):
        collection = self.db.filters
        doc = collection.find_one(dict(_id=ObjectId(filter_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.delete_one(dict(_id=ObjectId(filter_id)))

    def get_dataset_filter(self, email, dataset_id, filter_id):
        collection = self.db.filters
        doc = collection.find_one(dict(_id=ObjectId(filter_id)))
        self.get_dataset(email, doc['dataset_id'])
        return {'id': str(doc['_id']), 'dataset_id': doc['dataset_id'], 'name': doc['name'], 'value': doc['value'],
                'notes': doc['notes'], 'email': doc['email']}

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        self.get_dataset(email, dataset_id)
        collection = self.db.filters
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
        if filter_id is None:
            return str(collection.insert_one(entity_update).inserted_id)
        else:
            collection.update_one(dict(_id=ObjectId(filter_id)), {'$set': entity_update})
            return filter_id

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
            results.append(
                dict(id=str(doc['_id']), category=doc['category'], name=doc['name'], features=doc['features']))
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
        import datetime
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
                    submitted=doc.get('submitted')))
        return results

    def update_job(self, email, job_id, status, result):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)))
        self.get_dataset(email, doc['dataset_id'])
        collection.update_one(dict(_id=ObjectId(job_id)), {'$set': dict(status=status, result=result)})

    def delete_job(self, email, job_id):
        collection = self.db.jobs
        doc = collection.find_one(dict(_id=ObjectId(job_id)), dict(email=1))
        if doc['email'] == email:
            collection.delete_one(dict(_id=ObjectId(job_id)))
        else:
            raise InvalidUsage('Not authorized', 403)

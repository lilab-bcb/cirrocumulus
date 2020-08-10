import datetime
import json

from bson import ObjectId
from cirrocumulus.entity import Entity
from pymongo import MongoClient


class MongoDb:

    def __init__(self, db_uri, database='cirrocumulus', email=None):
        self.client = MongoClient(db_uri)
        self.db = self.client[database]
        self.email = email

    def server(self):
        d = {}
        if self.email is not None:
            d['email'] = self.email
        return d

    def category_names(self, dataset_id):
        collection = self.db.categories
        results = []
        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append({'category': doc['category'], 'dataset_id': doc['dataset_id'], 'original': doc['original'],
                            'new': doc['new']})
        return results

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        collection = self.db.categories
        key = str(dataset_id) + '-' + str(category) + '-' + str(original_name)

        if new_name == '':
            collection.delete_one(dict(_id=ObjectId(key)))
        else:
            collection.update_one(dict(_id=ObjectId(key)),
                {'$set': dict(category=category, dataset_id=dataset_id, original=original_name, new=new_name)},
                upsert=True)
            return key

    def user(self, email):
        collection = self.db.users
        collection.update_one(dict(_id=ObjectId(email)), {'$set': dict(last_login=datetime.datetime.now())},
            upsert=True)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        collection = self.db.datasets
        doc = collection.find_one(dict(_id=ObjectId(dataset_id)))
        if doc is None:
            raise ValueError('{} not found'.format(dataset_id))
        return Entity(str(doc['_id']), {'name': doc['name'],
                                        'url': doc['url'],
                                        'owner': 'owners' in doc and email in doc['owners']})

    def datasets(self, email):
        collection = self.db.datasets
        results = []
        for doc in collection.find():
            results.append({'id': str(doc['_id']), 'name': doc['name'],
                            'owner': 'owners' in doc and email in doc['owners']})
        return results

    def dataset_filters(self, email, dataset_id):
        collection = self.db.filters
        results = []

        for doc in collection.find(dict(dataset_id=dataset_id)):
            results.append(
                {'id': str(doc['_id']), 'dataset_id': doc['dataset_id'], 'name': doc['name'], 'value': doc['value'],
                 'notes': doc['notes'], 'email': doc['email']})
        return results

    def delete_dataset_filter(self, email, filter_id):
        collection = self.db.filters
        collection.delete_one(dict(_id=ObjectId(filter_id)))

    def get_dataset_filter(self, email, filter_id):
        collection = self.db.filters
        doc = collection.find_one(dict(_id=ObjectId(filter_id)))
        return {'id': str(doc['_id']), 'dataset_id': doc['dataset_id'], 'name': doc['name'], 'value': doc['value'],
                'notes': doc['notes'], 'email': doc['email']}

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
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
            collection.update_one(dict(_id=ObjectId(filter_id)), {'$set', entity_update})
            return filter_id

    def delete_dataset(self, email, dataset_id):
        collection = self.db.datasets
        collection.delete_one(dict(_id=ObjectId(dataset_id)))

    def upsert_dataset(self, email, dataset_id, dataset_name, url, readers):
        collection = self.db.datasets
        readers = set(readers)
        if email in readers:
            readers.remove(email)
        readers.add(email)
        update_dict = {'name': dataset_name,
                       'readers': list(readers),
                       'url': url}

        if dataset_id is None:  # new dataset
            update_dict['owners'] = [email]
            update_dict['_id'] = dataset_id
            return str(collection.insert_one(update_dict).inserted_id)
        else:
            collection.update_one(dict(_id=ObjectId(dataset_id)), {'$set', update_dict})
        return dataset_id

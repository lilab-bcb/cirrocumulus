import datetime
import json
import os

from cirrocumulus.abstract_db import AbstractDB
from cirrocumulus.envir import CIRRO_JOB_RESULTS, SERVER_CAPABILITY_EDIT_DATASET, SERVER_CAPABILITY_ADD_DATASET, \
    SERVER_CAPABILITY_DELETE_DATASET
from cirrocumulus.io_util import unique_id
from cirrocumulus.job_api import save_job_result_to_file
from cirrocumulus.util import get_fs


def create_dataset_meta(path):
    result = {'id': path, 'url': path, 'name': os.path.splitext(os.path.basename(path.rstrip('/')))[0],
              'description': ''}
    if os.path.basename(path).endswith('.json'):
        with open(path, 'rt') as f:
            result.update(json.load(f))

    return result


def write_json(json_data, json_path):
    if os.path.exists(os.path.dirname(os.path.abspath(json_path))):  # only support writing local files
        with open(json_path, 'wt') as f:
            json.dump(json_data, f)
    else:
        print('Skipping {}'.format(json_path))


class LocalDbAPI(AbstractDB):

    def __init__(self, paths):
        super().__init__()
        self.dataset_to_info = {}  # json_data, meta, json_path
        self.job_id_to_job = {}
        if os.environ.get(CIRRO_JOB_RESULTS) is not None:  # load saved on disk
            fs = get_fs(os.environ[CIRRO_JOB_RESULTS])
            fs.makedirs(os.environ[CIRRO_JOB_RESULTS], exist_ok=True)
            for url in fs.ls(os.environ[CIRRO_JOB_RESULTS]):
                if url.lower().endswith('.json.gz'):
                    import gzip
                    with gzip.open(fs.open(url)) as f:
                        d = json.load(f)
                        if 'id' in d:
                            d['url'] = url
                            self.job_id_to_job[d['id']] = d
                elif url.lower().endswith('.json'):
                    with fs.open(url) as f:
                        d = json.load(f)
                        if 'id' in d:
                            d['url'] = url
                            self.job_id_to_job[d['id']] = d

        for path in paths:
            json_data = {}
            basename = os.path.splitext(path)[0]
            old_path = basename + '_filters.json'
            json_path = basename + '.json'
            if os.path.exists(old_path) and os.path.getsize(old_path) > 0:
                with open(old_path, 'rt') as f:
                    json_data['filters'] = json.load(f)

            if os.path.exists(json_path) and os.path.getsize(json_path) > 0:
                with open(json_path, 'rt') as f:
                    try:
                        json_data.update(json.load(f))
                    except:
                        print('Unable to load {}'.format(json_path))
            meta = create_dataset_meta(path)
            if 'filters' not in json_data:
                json_data['filters'] = {}
            if 'views' not in json_data:
                json_data['views'] = {}
            if 'categories' not in json_data:
                json_data['categories'] = {}
            self.dataset_to_info[path] = dict(json_data=json_data, meta=meta, json_path=json_path)

    def capabilities(self):
        c = super().capabilities()
        c[SERVER_CAPABILITY_EDIT_DATASET] = False
        c[SERVER_CAPABILITY_ADD_DATASET] = False
        c[SERVER_CAPABILITY_DELETE_DATASET] = False
        return c

    def __delete_entity(self, dataset_id, entity_id, kind):
        json_data = self.dataset_to_info[dataset_id]['json_data']
        del json_data[kind][entity_id]
        write_json(json_data, self.dataset_to_info[dataset_id]['json_path'])

    def __get_entity(self, dataset_id, entity_id, kind):
        json_data = self.dataset_to_info[dataset_id]['json_data']
        return json_data[kind][entity_id]

    def __get_entity_list(self, dataset_id, kind):
        info = self.dataset_to_info.get(dataset_id)
        if info is None:
            return []
        json_data = info['json_data']
        results = []
        filters = json_data[kind]
        for key in filters:
            r = filters[key]
            r['id'] = key
            results.append(r)
        return results

    def __find_dataset_id(self, entity_id, kind):
        for dataset_id in self.dataset_to_info:
            info = self.dataset_to_info.get(dataset_id)
            if info is None:
                continue
            json_data = info['json_data']
            filters = json_data[kind]
            for key in filters:
                if key == entity_id:
                    return dataset_id
        raise ValueError('{} not found'.format(entity_id))

    def __upsert_entity(self, dataset_id, entity_id, kind, entity_dict):
        if entity_id is None:
            entity_id = unique_id()
        json_data = self.dataset_to_info[dataset_id]['json_data']
        entity = json_data[kind].get(entity_id)
        if entity is None:
            entity = {}
            json_data[kind][entity_id] = entity
        entity.update(entity_dict)
        write_json(json_data, self.dataset_to_info[dataset_id]['json_path'])
        return entity_id

    def user(self, email):
        return dict()

    def datasets(self, email):
        results = []
        for key in self.dataset_to_info:
            results.append(self.dataset_to_info[key]['meta'])
        return results

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        info = self.dataset_to_info.get(dataset_id)
        if info is None:  # on the fly
            result = {'id': dataset_id, 'url': dataset_id, 'name': os.path.splitext(os.path.basename(dataset_id))[0],
                      'description': ''}
        else:
            result = info['meta']
            result['id'] = dataset_id

        return result

    def category_names(self, email, dataset_id):
        info = self.dataset_to_info.get(dataset_id)
        if info is None:
            return {}
        json_data = info['json_data']
        categories = json_data['categories']
        results = {}
        for category_name in categories:
            category = categories[category_name]
            cat_results = {}
            results[category_name] = cat_results
            for original_value in category:
                val = category[original_value]
                cat_results[original_value] = dict(originalValue=original_value, newValue=val.get('new'),
                                                   positiveMarkers=val.get('positiveMarkers'),
                                                   negativeMarkers=val.get('negativeMarkers'),
                                                   color=val.get('color'))
        return results

    def upsert_category_name(self, email, dataset_id, category, original_value, update):
        json_data = self.dataset_to_info[dataset_id]['json_data']
        category_entity = json_data['categories'].get(category)
        if category_entity is None:
            category_entity = {}
            json_data['categories'][category] = category_entity
        if 'newValue' in update:
            update['new'] = update.pop('newValue')
        entity = category_entity.get(original_value)
        if entity is None:
            entity = {}
        entity.update(update)
        category_entity[original_value] = entity
        if dataset_id is not None:
            entity['dataset_id'] = dataset_id
        if email is not None:
            entity['email'] = email
        write_json(json_data, self.dataset_to_info[dataset_id]['json_path'])

    def get_feature_sets(self, email, dataset_id):
        info = self.dataset_to_info.get(dataset_id)
        if info is None:
            return []
        json_data = info['json_data']
        return json_data.get('markers', [])

    def delete_feature_set(self, email, dataset_id, set_id):
        json_data = self.dataset_to_info[dataset_id]['json_data']
        markers = json_data['markers']
        for i in range(len(markers)):
            if markers[i]['id'] == set_id:
                markers.pop(i)
                break
        write_json(json_data, self.dataset_to_info[dataset_id]['json_path'])

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        if set_id is None:
            set_id = unique_id()
        else:
            self.delete_feature_set(email=email, dataset_id=dataset_id, set_id=set_id)
        json_data = self.dataset_to_info[dataset_id]['json_data']
        markers = json_data.get('markers')
        if markers is None:
            markers = []
            json_data['markers'] = markers
        markers.append(dict(id=set_id, features=features, name=name, category=category))
        write_json(json_data, self.dataset_to_info[dataset_id]['json_path'])
        return set_id

    def dataset_views(self, email, dataset_id):
        return self.__get_entity_list(dataset_id=dataset_id, kind='views')

    def delete_dataset_view(self, email, view_id):
        dataset_id = self.__find_dataset_id(view_id, 'views')
        return self.__delete_entity(dataset_id=dataset_id, entity_id=view_id, kind='views')

    def get_dataset_view(self, email, view_id):
        dataset_id = self.__find_dataset_id(view_id, 'views')
        return self.__get_entity(dataset_id=dataset_id, entity_id=view_id, kind='views')

    def upsert_dataset_view(self, email, dataset_id, view):
        view['last_updated'] = str(datetime.datetime.utcnow())
        if 'value' in view:
            view['value'] = json.dumps(view['value'])
        if email is not None:
            view['email'] = email
        if dataset_id is not None:
            view['dataset_id'] = dataset_id
        view_id = view.pop('id') if 'id' in view else None
        view_id = self.__upsert_entity(dataset_id=dataset_id, entity_id=view_id, kind='views', entity_dict=view)
        return dict(id=view_id, last_updated=view['last_updated'])

    def create_job(self, email, dataset_id, job_name, job_type, params):
        job_id = unique_id()
        self.job_id_to_job[job_id] = dict(id=job_id, dataset_id=dataset_id, name=job_name, type=job_type, params=params,
                                          status=None, submitted=str(datetime.datetime.utcnow()))
        return job_id

    def get_job(self, email, job_id, return_type):
        job = self.job_id_to_job.get(job_id)
        if job is None:
            return None
        if return_type == 'result':
            return job['result'] if 'result' in job else job
        elif return_type == 'status':
            return dict(status=job['status'])
        elif return_type == 'params':
            return dict(params=job.get('params'))

    def get_jobs(self, email, dataset_id):
        results = []
        for job in self.job_id_to_job.values():
            if dataset_id == job['dataset_id']:
                results.append(dict(id=job['id'], name=job['name'], type=job['type'], status=job['status'],
                                    submitted=job['submitted']))
        return results

    def delete_job(self, email, job_id):
        job = self.job_id_to_job.pop(job_id)
        if 'url' in job and os.path.exists(job['url']):
            os.remove(job['url'])

    def update_job(self, email, job_id, status, result):
        job = self.job_id_to_job[job_id]
        job['status'] = status
        if result is not None:
            if os.environ.get(CIRRO_JOB_RESULTS) is not None:  # save to directory
                result.update(job)
                save_job_result_to_file(result, job_id)
                job['result'] = result
            else:
                from cirrocumulus.util import to_json
                job['result'] = to_json(result)

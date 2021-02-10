def get_email_domain(email):
    at_index = email.find('@')
    domain = None
    if at_index != -1:
        domain = email[at_index + 1:]
    return domain


def load_dataset_schema(url):
    import os
    from urllib.parse import urlparse
    import fsspec
    import json

    def get_extension(path):
        name, ext = os.path.splitext(path)
        if ext == '.gz':
            name, ext = os.path.splitext(name)
            if ext == '.json':
                ext = '.json.gz'
        return ext

    extension = get_extension(url)
    json_schema = None
    if extension in ['.json', '.json.gz', '']:
        pr = urlparse(url)
        fs = fsspec.filesystem(pr.scheme if not pr.scheme == '' else 'file')
        if extension == '':
            url = os.path.join(url, 'index.json.gz')
            extension = get_extension(url)
        if extension == '.json.gz':
            import gzip
            with gzip.open(fs.open(url)) as f:
                json_schema = json.load(f)
        else:
            with fs.open(url) as f:
                json_schema = json.load(f)
    return json_schema


class DatabaseAPI:

    def __init__(self):
        self.provider = None

    def server(self):
        return self.provider.server()

    def user(self, email):
        return self.provider.user(email)

    def datasets(self, email):
        return self.provider.datasets(email)

    def category_names(self, email, dataset_id):
        return self.provider.category_names(email, dataset_id)

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        return self.provider.upsert_category_name(email, category, dataset_id, original_name, new_name)

    def delete_dataset(self, email, dataset_id):
        return self.provider.delete_dataset(email, dataset_id)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        return self.provider.get_dataset(email, dataset_id, ensure_owner)

    def upsert_dataset(self, email, dataset_id, dataset_name, url, readers, description, title, species):
        return self.provider.upsert_dataset(email=email, dataset_id=dataset_id, dataset_name=dataset_name, url=url,
            readers=readers, description=description, title=title, species=species)

    def dataset_filters(self, email, dataset_id):
        return self.provider.dataset_filters(email, dataset_id)

    def get_dataset_filter(self, email, dataset_id, filter_id):
        return self.provider.get_dataset_filter(email, dataset_id, filter_id)

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        return self.provider.upsert_dataset_filter(email, dataset_id, filter_id, filter_name, filter_notes,
            dataset_filter)

    def delete_dataset_filter(self, email, dataset_id, filter_id):
        return self.provider.delete_dataset_filter(email, dataset_id, filter_id)

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        return self.provider.upsert_feature_set(email, dataset_id, set_id, category, name, features)

    def get_feature_sets(self, email, dataset_id):
        return self.provider.get_feature_sets(email, dataset_id)

    def delete_feature_set(self, email, dataset_id, set_id):
        return self.provider.delete_feature_set(email, dataset_id, set_id)

    def create_job(self, email, dataset_id, job_name, job_type, params):
        return self.provider.create_job(email, dataset_id, job_name, job_type, params)

    def update_job(self, email, job_id, status, result):
        return self.provider.update_job(email, job_id, status, result)

    def get_job(self, email, job_id, return_result):
        return self.provider.get_job(email, job_id, return_result)

    def get_jobs(self, email, dataset_id):
        return self.provider.get_jobs(email, dataset_id)

    def delete_job(self, email, job_id):
        return self.provider.delete_job(email, job_id)

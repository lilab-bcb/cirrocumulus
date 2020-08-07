class DatabaseAPI:

    def __init__(self):
        self.provider = None

    def server(self):
        return self.provider.server()

    def user(self, email):
        return self.provider.user(email)

    def datasets(self, email):
        return self.provider.datasets(email)

    def category_names(self, dataset_id):
        return self.provider.category_names(dataset_id)

    def upsert_category_name(self, email, category, dataset_id, original_name, new_name):
        return self.provider.upsert_category_name(email, category, dataset_id, original_name, new_name)

    def delete_dataset(self, email, dataset_id):
        return self.provider.delete_dataset(email, dataset_id)

    def get_dataset(self, email, dataset_id, ensure_owner=False):
        return self.provider.get_dataset(email, dataset_id, ensure_owner)

    def upsert_dataset(self, email, dataset_id, dataset_name, url, readers):
        return self.provider.upsert_dataset(email, dataset_id, dataset_name, url, readers)

    def dataset_filters(self, email, dataset_id):
        return self.provider.dataset_filters(email, dataset_id)

    def get_dataset_filter(self, email, filter_id):
        return self.provider.get_dataset_filter(email, filter_id)

    def upsert_dataset_filter(self, email, dataset_id, filter_id, filter_name, filter_notes, dataset_filter):
        return self.provider.upsert_dataset_filter(email, dataset_id, filter_id, filter_name, filter_notes,
            dataset_filter)

    def delete_dataset_filter(self, email, filter_id):
        return self.provider.delete_dataset_filter(email, filter_id)

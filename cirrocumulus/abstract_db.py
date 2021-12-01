# Abstract database class supporting datasets, categories, views, feature sets, users, jobs, and job results
# Only methods server and datasets must be supported in 'client' mode
import os

from cirrocumulus.envir import SERVER_CAPABILITY_RENAME_CATEGORIES, SERVER_CAPABILITY_JOBS, \
    SERVER_CAPABILITY_FEATURE_SETS, SERVER_CAPABILITY_LINKS, SERVER_CAPABILITY_EDIT_DATASET, \
    SERVER_CAPABILITY_ADD_DATASET, SERVER_CAPABILITY_DELETE_DATASET


def to_bool(s):
    return s.lower() in ['true', '1']


class AbstractDB:

    def __init__(self):
        """Initializes the object

        """

    def capabilities(self):  # allow everything
        c = {}
        c[SERVER_CAPABILITY_RENAME_CATEGORIES] = to_bool(os.environ.get(SERVER_CAPABILITY_RENAME_CATEGORIES, 'True'))
        c[SERVER_CAPABILITY_JOBS] = to_bool(os.environ.get(SERVER_CAPABILITY_JOBS, 'True'))
        c[SERVER_CAPABILITY_FEATURE_SETS] = to_bool(os.environ.get(SERVER_CAPABILITY_FEATURE_SETS, 'True'))
        c[SERVER_CAPABILITY_LINKS] = to_bool(os.environ.get(SERVER_CAPABILITY_LINKS, 'True'))
        c[SERVER_CAPABILITY_EDIT_DATASET] = to_bool(os.environ.get(SERVER_CAPABILITY_EDIT_DATASET, 'True'))
        c[SERVER_CAPABILITY_ADD_DATASET] = to_bool(os.environ.get(SERVER_CAPABILITY_ADD_DATASET, 'True'))
        c[SERVER_CAPABILITY_DELETE_DATASET] = to_bool(os.environ.get(SERVER_CAPABILITY_DELETE_DATASET, 'True'))
        return c

    def datasets(self, email):
        """ Gets list of available datasets

        Args:
            email: User email or None

        Returns:
            A list of dicts. Example:
            [{"id": "dataset_id",
             "name": "dataset_name"
             "title": "dataset_title",
             "owner": false,
             "url": "gs://xx/my_dataset".
             "species": "Mus Musculus"
            }]
        """
        raise NotImplementedError()

    def category_names(self, email, dataset_id):
        """Gets a dict of renamed categories for a dataset

        Args:
              email: User email or None
              dataset_id: Dataset id

        Returns:
            A dict that maps category->category_value->{}. Example:

            {"louvain":{"1":{"color":"red"}}
        """
        raise NotImplementedError()

    def upsert_category_name(self, email, dataset_id, category, original_value, update):
        """ Upserts a category name.

        Args:
             email: User email or None
             category: Category in dataset (e.g. louvain)
             dataset_id: Dataset id
             original_value: Original category value (e.g. "1")
             update: Update object optionally containing color, positiveMarkers, negativeMarkers, newValue
        """
        raise NotImplementedError()

    def user(self, email):
        """ Gets metadata about a user

        Args:
           email: User email

        Returns:
            A dict with the keys id and importer. An importer can add datasets to cirrocumulus. Example:
            {"id":"user_id",
            "importer":false
            }
        """
        raise NotImplementedError()

    # views
    def dataset_views(self, email, dataset_id):
        """ Gets list of saved dataset views (saved visualization states)

       Args:
             email: User email or None
             dataset_id: Dataset id

        Returns:
            List of dicts. Example:
            [{"id": "view id",
              "name": "view name"}]
       """
        raise NotImplementedError()

    def delete_dataset_view(self, email, view_id):
        """ Delete a saved view

        Args:
            email: User email or None
            view_id: View id

        """
        raise NotImplementedError()

    def get_dataset_view(self, email, view_id):
        """ Gets detailed information for a saved dataset view

        Args:
            email: User email or None
            view_id: View id

        Returns:
            List of dicts containing id and name. Example:
            [{"id": "view id",
              "name": "view name",
              "value": "JSON encoded state"
              "notes": "view notes"
              "email": "View creator email"
        """

        raise NotImplementedError()

    def upsert_dataset_view(self, email, dataset_id, view):
        """ Upserts a dataset view
        View should have id (for update), name, value, and any other additional fields to store

        Args:
              email: User email or None
              dataset_id: Dataset id
              view: View to upsert

         Returns:
            dict with id, last_updated
        """
        raise NotImplementedError()

    def delete_dataset(self, email, dataset_id):
        """ Deletes a dataset

        Args:
            email: User email or None
            dataset_id: Dataset id
        """
        raise NotImplementedError()

    def upsert_dataset(self, email, readers, dataset):
        """ Upserts a dataset. If dataset.id is None then a new dataset is inserted.
        Dataset should have name, url, description, title, species, and any other additional fields to store
        Args:
              email: User email or None
              readers: List of allowed readers
              dataset: Dataset to upsert

         Returns:
            Upserted dataset id
        """
        raise NotImplementedError()

    def get_feature_sets(self, email, dataset_id):
        """ Gets saved feature sets

        Args:
              email: User email or None
              dataset_id: Dataset id

        Returns:
            List of dicts. Example:
            [{"id": "set_id",
              "category": "set category",
              "name": "set name",
              "features: ["gene1", "gene2"]}]
        """
        raise NotImplementedError()

    def delete_feature_set(self, email, dataset_id, set_id):
        """ Deletes a saved feature set

        Args:
              email: User email or None
              dataset_id: Dataset id
              set_id: Feature set id
        """
        raise NotImplementedError()

    def upsert_feature_set(self, email, dataset_id, set_id, category, name, features):
        """ Upserts a feature set

        Args:
            email: User email or None
            dataset_id: Dataset id
            set_id: Feature set id
            category: Set category
            name: Set name
            features: List of features

        Returns:
            Upserted id
        """
        raise NotImplementedError()

    def create_job(self, email, dataset_id, job_name, job_type, params):
        """ Creates a job

        Args:
         email: User email or None
         dataset_id: Dataset id
         job_name: Job name
         job_type: Job type
         params: JSON encoded job params

        Returns:
         job id
      """
        raise NotImplementedError()

    def get_job(self, email, job_id, return_type):
        """ Gets a job

       Args:
          email: User email or None
          job_id: Job id
          return_type: One of "result", "status", or "params"

       Returns:
          The job
      """
        raise NotImplementedError()

    def get_jobs(self, email, dataset_id):
        """ Gets a list of all jobs for a dataset.

        Args:
            email: User email or None
            dataset_id: Dataset id

        Returns:
            List of jobs
      """
        raise NotImplementedError()

    def delete_job(self, email, job_id):
        """ Deletes a job.

      Args:
          email: User email or None
          job_id: Job id
      """
        raise NotImplementedError()

    def update_job(self, email, job_id, status, result):
        """ Updates job info.

          Args:
              email: User email or None
              job_id: Job id
              status: Job status
              result: Job result
          """
        raise NotImplementedError()

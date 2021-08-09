# Abstract database class supporting datasets, categories, views, feature sets, users, jobs, and job results
# Only methods server and datasets must be supported in 'client' mode
from cirrocumulus.envir import *
from cirrocumulus.io_util import unique_id


class AbstractDB:

    def __init__(self):
        """Initializes the object

        """
        self.job_id_to_job = {}

    def capabilities(self):  # allow everything
        c = {}
        c[SERVER_CAPABILITY_RENAME_CATEGORIES] = True
        c[SERVER_CAPABILITY_JOBS] = True
        c[SERVER_CAPABILITY_SAVE_FEATURE_SETS] = True
        c[SERVER_CAPABILITY_SAVE_LINKS] = True
        c[SERVER_CAPABILITY_EDIT_DATASET] = True
        c[SERVER_CAPABILITY_ADD_DATASET] = True
        c[SERVER_CAPABILITY_DELETE_DATASET] = True
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
        """Gets a list of renamed category names for a dataset

        Args:
              email: User email or None
              dataset_id: Dataset id

        Returns:
            A list of dicts representing renamed categories. Example:

            {"category":"louvain",
            "dataset_id":"1",
            "original":"1",
            "new":"my cell type"}
        """
        raise NotImplementedError()

    def upsert_category_name(self, email, category, dataset_id, original_name, prior_name, new_name):
        """ Upserts a category name.

        Args:
             email: User email or None
             category: Category in dataset (e.g. louvain)
             dataset_id: Dataset id
             original_name: Original category name (e.g. "1")
             prior_name: Prior name (e.g. "cell type a")
             new_name: New name (e.g. "cell type b")
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

    def delete_dataset_view(self, email, dataset_id, view_id):
        """ Delete a saved view

        Args:
            email: User email or None
            dataset_id: Dataset id
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

    def upsert_dataset_view(self, email, dataset_id, view_id, name, value):
        """ Upserts a dataset view

        Args:
              email: User email or None
              dataset_id: Dataset id
              view_id: View id or None to create new view
              name: View name
              value: JSON encoded state

         Returns:
            Upserted view id
        """
        raise NotImplementedError()

    def delete_dataset(self, email, dataset_id):
        """ Deletes a dataset

        Args:
            email: User email or None
            dataset_id: Dataset id
        """
        raise NotImplementedError()

    def upsert_dataset(self, email, dataset_id, dataset_name=None, url=None, readers=None, description=None, title=None,
                       species=None):
        """ Upserts a dataset

        Args:
              email: User email or None
              dataset_id: Dataset id
              dataset_name: Name
              url: URL
              readers: List of allowed readers
              description: Description
              title: Title
              species: Species

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
        import datetime
        job_id = unique_id()
        self.job_id_to_job[job_id] = dict(id=job_id, dataset_id=dataset_id, name=job_name, type=job_type, params=params,
                                          status=None, result=None, submitted=datetime.datetime.utcnow())
        return job_id

    def get_job(self, email, job_id, return_result):
        """ Gets a job

       Args:
          email: User email or None
          job_id: Job id
          return_result: Whether to return the job result or status only

       Returns:
          The job
      """
        job = self.job_id_to_job[job_id]
        if return_result:
            return job['result']
        return dict(id=job['id'], name=job['name'], type=job['type'], status=job['status'], submitted=job['submitted'])

    def get_jobs(self, email, dataset_id):
        """ Gets a list of all jobs for a dataset.

        Args:
            email: User email or None
            dataset_id: Dataset id

        Returns:
            List of jobs
      """
        results = []
        for job in self.job_id_to_job.values():
            results.append(dict(id=job['id'], name=job['name'], type=job['type'], status=job['status'],
                                submitted=job['submitted']))
        return results

    def delete_job(self, email, job_id):
        """ Deletes a job.

      Args:
          email: User email or None
          job_id: Job id
      """
        del self.job_id_to_job[job_id]

    def update_job(self, email, job_id, status, result):
        """ Updates job info.

          Args:
              email: User email or None
              job_id: Job id
              status: Job status
              result: Job result
          """
        job = self.job_id_to_job[job_id]
        job['status'] = status
        if result is not None:
            from cirrocumulus.util import to_json
            job['result'] = to_json(result)

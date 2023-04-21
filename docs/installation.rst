Installation
-------------

Standalone Mode
^^^^^^^^^^^^^^^^

- Install using pip::

    pip install cirrocumulus

- Launch via the command line::

    cirro launch <path_to_dataset>

- Full list of command line options:

    .. argparse::
       :ref: cirrocumulus.launch.create_parser
       :prog: cirro launch


Server Mode
^^^^^^^^^^^^^^

Cirrocumulus can also be run in `server` mode in order to serve multiple users and datasets securely.
The cirrocumulus server can be deployed on a cloud VM, an on-premise machine, or on Google App Engine.

- Install cirrocumulus using pip or docker

- Optional additional setup to enable authentication and authorization via Okta or Google:

    - Install additional libraries:

        - Google: `pip install google-auth`
        - Okta: `pip install okta-jwt-verifier`

    - Set environment variables:

        - CIRRO_AUTH_CLIENT_ID: to your Okta or Google client id
        - CIRRO_AUTH_PROVIDER: to either okta or Google.
        - CIRRO_AUTH_ISSUER (okta only). The URL of the authorization server that will perform authentication.
          All Developer Accounts have a "default" authorization server.
          The issuer is a combination of your Org URL (found in the upper right of the console home page)
          and /oauth2/default. For example, https://dev-1234.oktapreview.com/oauth2/default.

    - See `Okta documentation`_ for creating custom app integrations with Okta.

    - Visit the `Google OAuth 2.0 documentation`_ to obtain OAuth 2.0 credentials for Google.

    - Please note that https is required if using Okta or Google authentication

- Additional libraries needed for cloud storage:

    - Amazon S3: `pip install s3fs`
    - Google Cloud Storage: `pip install gcsfs`
    - Microsoft Azure: `pip install adlfs`

- Install MongoDB_ and start the MongoDB server. Note that MongoDB compatible databases such as DocumentDB_ can also be used.

- Start the server via the command line::

    cirro serve


- Use the :ref:`prepare_data<prepare_data>` command to freeze an h5ad, loom, or Seurat file in cirrocumulus format.

- Add a dataset and optionally share with dataset with collaborators. If you enabled authentication, then no users are allowed to add datasets to cirrocumulus.
  Set the property "importer" to true on an entry in the users collection to enable that user to import datasets. For example, the following screenshot in `MongoDB Compass`_ shows that the user with the email address `me@gmail.com`, is allowed to add datasets to cirrocumulus:

    .. image:: images/mongodb.png


- You can programmatically add a dataset by posting to the /api/dataset endpoint::

    curl http://localhost:5000/api/dataset -X POST -F 'name=my_name' -F 'url=data/my_dataset_path' -F 'description=my_desc'  -F 'species=Mus musculus'

- Additional customization via environment variables:

    - CIRRO_MOUNT: For mounting a bucket locally. Comma separated string of bucket:local_path. Example s3://foo/bar:/fsx
    - CIRRO_SPECIES: Path to JSON file for species list when adding new dataset
    - CIRRO_MIXPANEL: Mixpanel_ project token for event tracking. Currently, only the open dataset event is supported.

- Optionally, set the default view for a dataset by adding the field "defaultView" to your dataset entry in the database.
  You can configure the cirrocumulus state in the app, then use "Copy Link" to get the JSON configuration for "defaultView". Example:

    .. image:: images/default_view.png

- Full list of command line options:

    .. argparse::
       :ref: cirrocumulus.serve.create_parser
       :prog: cirro serve

Static Website
^^^^^^^^^^^^^^^^

- Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus.git

- Change to cirrocumulus directory::

    cd cirrocumulus

- Install typescript::

    yarn global add typescript

- Install JavaScript dependencies::

    yarn install

- Prepare dataset(s) in jsonl format::

    cirro prepare_data pbmc3k.h5ad --format jsonl

- Build JavaScript::

    REACT_APP_STATIC=true yarn build

- Create the file datasets.json in the build directory::


    [
        {
            "id": "pbmc3k",
            "name": "pbmc3k",
            "url": "pbmc3k/pbmc3k.jsonl"
        }
    ]


- Move your dataset files to build::

    mv pbmc3k build

- Test locally::

    cd build ; npx http-server .

- Host the build directory on your static website hosting service (e.g. `Amazon S3`_, `Google Cloud Storage`_)

.. _prepare_data:

Prepare Data
^^^^^^^^^^^^^^
The `prepare_data` command is used to freeze an h5ad, loom, or Seurat (RDS) file in cirrocumulus format.
The cirrocumulus format allows efficient partial dataset retrieval over a network (e.g Google bucket) using limited memory.

- Example::

    cirro prepare_data pbmc3k.h5ad


- Full list of command line options:

    .. argparse::
       :ref: cirrocumulus.prepare_data.create_parser
       :prog: cirro prepare_data


Developer Instructions
^^^^^^^^^^^^^^^^^^^^^^^^

- Create a new conda environment::

    conda create --name cirrocumulus-dev

- Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus.git

- Change to cirrocumulus directory::

    cd cirrocumulus


- Install::

    pip install --upgrade pip
    pip install -e .[dev,test]
    pre-commit install
    yarn global add typescript
    yarn install
    yarn build
    pip install -e .

- Install additional optional Python dependencies::

    pip install s3fs

- Create an example h5ad file in ./data/pbmc3k_processed.h5ad::

    import scanpy as sc
    sc.datasets.pbmc3k_processed()

- Launch cirrocumulus with the --no-open flag::

    cirro launch ./data/pbmc3k_processed.h5ad --no-open

- Alternatively, launch the cirrocumulus server (see :ref:`prepare_data<prepare_data>`)::

    cirro serve

- Run JavaScript server in development mode::

    yarn start

- Navigate to http://localhost:3000

- In order to run End to End tests (yarn e2e), please install GraphicsMagick (brew install graphicsmagick on Mac)

- Testing::

    yarn e2e
    yarn test
    pytest

- Build JavaScript front-end for deployment::

    yarn build



.. _app.yaml: https://cloud.google.com/appengine/docs/standard/python3/config/appref
.. _Google Cloud SDK: https://cloud.google.com/sdk/install
.. _App Engine: https://cloud.google.com/appengine/docs/
.. _Node.js: https://nodejs.org/
.. _ngrok: https://ngrok.com/
.. _Terra: https://app.terra.bio/
.. _MongoDB: https://www.mongodb.com/
.. _Google API Console: https://console.developers.google.com/
.. _gcsfuse: https://github.com/GoogleCloudPlatform/gcsfuse/
.. _MongoDB Compass: https://www.mongodb.com/products/compass
.. _Amazon S3: https://docs.aws.amazon.com/AmazonS3/latest/userguide/WebsiteHosting.html
.. _Google Cloud Storage: https://cloud.google.com/storage/docs/hosting-static-website-http
.. _Mixpanel: https://mixpanel.com/
.. _Okta documentation: https://help.okta.com/en/prod/Content/Topics/Apps/Apps_App_Integration_Wizard.htm
.. _Google OAuth 2.0 documentation: https://support.google.com/cloud/answer/6158849
.. _Cloud Build API: https://console.cloud.google.com/flows/enableapi?apiid=cloudbuild.googleapis.com
.. _DocumentDB: https://aws.amazon.com/documentdb/

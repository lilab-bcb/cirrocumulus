Installation
-------------

pip
^^^^^

- Install using pip::

    pip install cirrocumulus

- Launch via the command line::

    cirro launch <path_to_dataset>


Terra_ notebook
^^^^^^^^^^^^^^^^
- Click ``Open Terminal`` to connect to your running VM
- Install cirrocumulus via pip if it was not installed in your docker image
- Download your dataset to your running VM using gsutil::

    gsutil -m cp gs://fc-000/test.h5ad .

- Launch cirrocumulus via the command line in the background::

    cirro launch test.h5ad &

- Install ngrok_::

    wget https://bin.equinox.io/c/4VmDzA7iaHb/ngrok-stable-linux-amd64.zip \
    && unzip ngrok-stable-linux-amd64.zip \
    && rm -f ngrok-stable-linux-amd64.zip

- Use ngrok_ to expose cirrocumulus publicly::

    ./ngrok http 5000

After you start ngrok, it will display a UI in your terminal with the public URL of your tunnel.

- Navigate to your public URL in your browser.

Docker
^^^^^^^^

- Launch using docker::

    docker run -it -p 5000:5000 --rm -v `pwd`:/data cumulusprod/cirrocumulus cirro launch /data/dataset1.h5ad --host 0.0.0.0




Build From Source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus


- Install JavaScript dependencies, build the client, and install cirrocumulus Python module in editable mode::

    ./build.sh

- Launch cirrocumulus via the command line::

    cirro launch test.h5ad

Server Mode
^^^^^^^^^^^^^^

Cirrocumulus can also be run in `server` mode in order to serve multiple users and datasets.
The server can be deployed on a cloud VM, an on-premise machine, or on Google App Engine.


Prepare Data
^^^^^^^^^^^^^^
The command `cirro prepare_data` can be used to prepare a cirrocumulus formatted dataset from a h5ad, loom, or Seurat file.
The cirrocumulus format allows efficient fetching of portions of a dataset over a network (e.g. Google or S3 bucket). Example::

    cirro prepare_data pbmc3k_final.rds


Google App Engine
^^^^^^^^^^^^^^^^^^^


- Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus

- Change your current working directory to cirrocumulus::

    cd cirrocumulus

- Install `Node.js`_.

- Install TypeScript::

    npm install -g typescript

- Install JavaScript dependencies::

    npm i


- Build the client::

    npm run-script build

- Create or use an existing GCP project

- Create App Engine by navigating to App Engine > Dashboard. You may choose the region where your application is hosted.
  Select the Python 3 Standard Environment.

- Obtain OAuth 2.0 credentials.
    - Create an OAuth client id. Set the OAuth consent screen application name and add your server URL to the list of “Authorized domains”. Your server URL is \https://<PROJECT>.appspot.com.
    - Go to Credentials and click “Create Credentials > OAuth client ID”. Enter “Web application” for “Application Type”
      and your server URL for “Authorized JavaScript origins”. Click “Create” to create the credentials.

- Replace CIRRO_AUTH_CLIENT_ID in app.yaml with your OAuth client id.

- Install the `Google Cloud SDK`_. Type ``gcloud init`` in your terminal if this is your first time using the Google Cloud SDK.

- Deploy the application using the command below. Remember to replace <PROJECT> with your project ID. Your project is available at \https://<PROJECT>.appspot.com.::

    gcloud app deploy app.yaml --project=<PROJECT>

- Go to \https://<PROJECT>.appspot.com in your web browser and login.

    - By default, no one is allowed to add datasets to your application.
    - In Google Console, navigate to Data Store > Entities and click on your email address. Add the property ``importer`` of type ``boolean`` and set it to ``true``.
    - Go back to \https://<PROJECT>.appspot.com and start adding datasets.

- Read more about App Engine in the `App Engine`_ documentation.


Server
^^^^^^^^

- Install cirrocumulus using pip or docker

- Optionally visit the `Google API Console`_ to obtain OAuth 2.0 credentials.

    - Create an OAuth client id. Set the OAuth consent screen application name and add your server URL to the list of “Authorized domains”
    - Go to Credentials and click “Create Credentials > OAuth client ID”. Enter “Web application” for “Application Type”
      and your server URL for “Authorized JavaScript origins”. Click “Create” to create the credentials.

- Install MongoDB_ and start the MongoDB server

- Start the server via the command line::

    cirro serve

- Import a dataset

.. _Google Cloud SDK: https://cloud.google.com/sdk/install
.. _App Engine: https://cloud.google.com/appengine/docs/
.. _Node.js: https://nodejs.org/
.. _ngrok: https://ngrok.com/
.. _Terra: https://app.terra.bio/
.. _MongoDB: https://www.mongodb.com/
.. _Google API Console: https://console.developers.google.com/
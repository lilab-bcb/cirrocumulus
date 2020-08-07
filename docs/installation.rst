Installation
-------------

pip
^^^^^

- Install using pip::

    pip install cirrocumulus

-  Launch via the command line::

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


Google Cloud Platform (GCP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus

-  Change your current working directory to cirrocumulus::

    cd cirrocumulus

-  Install `Node.js`_.

-  Install TypeScript::

    npm install -g typescript

-  Install JavaScript dependencies::

    npm i


-  Build the client::

    npm run-script build

-  Create or use an existing GCP project

-  Create an OAuth client id

   -  In Google Console, navigate to APIs and Services > OAuth consent
      screen. Set the OAuth consent screen application name and add
      <PROJECT>.appspot.com to the list of “Authorized domains”
   -  Go to Credentials and click “Create Credentials > OAuth client
      ID”. Enter “Web application” for “Application Type” and
      https://<PROJECT>.appspot.com for “Authorized JavaScript origins”.
      Click “Create” to create the credentials.
   -  Copy OAuth client id into ``cirrocumulus/cirrocumulus/config.json``.

-  Create App Engine by navigating to App Engine > Dashboard. You may
   choose the region where your application is hosted. Select the Python
   3 Standard Environment.
-  Install the `Google Cloud SDK`_. Type ``gcloud init`` in your terminal if this is your
   first time using the Google Cloud SDK.
-  Deploy the application using the command below. Remember to replace
   <PROJECT> with your project ID.::

    gcloud app deploy app.yaml --project=<PROJECT>

   Your project is available at https://<PROJECT>.appspot.com.

-  Go to https://<PROJECT>.appspot.com in your web browser and login.

   -  By default, no one is allowed to add datasets to your application.
   -  In Google Console, navigate to Data Store > Entities and click on
      your email address. Add the property ``importer`` of type ``boolean``
      and set it to ``true``.
   -  Go back to https://<PROJECT>.appspot.com and start adding datasets.

-  Read more about `App Engine`_ in the `App Engine`_ documentation.


Build From Source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus


-  Install JavaScript dependencies, build the client, and install cirrocumulus Python module in editable mode::

    ./build.sh

- Launch cirrocumulus via the command line::

    cirro launch test.h5ad


.. _Google Cloud SDK: https://cloud.google.com/sdk/install
.. _App Engine: https://cloud.google.com/appengine/docs/
.. _Node.js: https://nodejs.org/
.. _ngrok: https://ngrok.com/
.. _Terra: https://app.terra.bio/
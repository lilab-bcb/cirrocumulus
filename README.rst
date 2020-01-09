cirrocumulus
===============
Interactive single cell visualization application.
 - Run on Google App Engine or your desktop
 - Scales to > 6 million cells on commodity hardware
 - Share datasets with collaborators
 - Interactively create and share cell annotations
 - Explore multiple features and embeddings simultaneously

Local Usage
-----------

- Install using pip::

   pip install cirrocumulus

-  Launch via the command line::

    cirro launch <path to dataset>

Google Cloud Platform (GCP) Deployment Instructions
-----------------------------------------------------

-  Clone the cirrocumulus repository::

    git clone https://github.com/klarman-cell-observatory/cirrocumulus

-  Change your current working directory to cirrocumulus::

    cd cirrocumulus

-  Install `Node.js`_.

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
   -  Copy OAuth client id into ``cirrocumulus/cirro/config.json``.

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

-  Read more about `App Engine`_, such as how you can limit spending.

Developer Instructions
--------------------------

-  Install JavaScript dependencies::

    npm i

-  Build the client::

    npm run-script build

-  Install Python module in editable mode::

    pip install -r requirements.txt -e .

-  Add http://localhost:5000 to your Web application Outh client ID
   authorized JavaScript origins at APIs and Services > Credentials
-  Download the App Engine service account JSON key from IAM & admin > Service accounts (DO NOT SHARE THIS!)
   and set the environment variable GOOGLE_APPLICATION_CREDENTIALS::

    export GOOGLE_APPLICATION_CREDENTIALS=“/home/user/Downloads/service-account-file.json”


.. image:: https://travis-ci.org/klarman-cell-observatory/cirrocumulus.svg?branch=master
    :target: https://travis-ci.org/klarman-cell-observatory/cirrocumulus

.. _Google Cloud SDK: https://cloud.google.com/sdk/install
.. _App Engine: https://cloud.google.com/appengine/docs/
.. _Node.js: https://nodejs.org/
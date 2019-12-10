import os

from flask import Flask, send_from_directory

from cirro.api import blueprint, auth_api, database_api, dataset_api
from cirro.firestore_datastore import FirestoreDatastore
from cirro.google_auth import GoogleAuth
from cirro.parquet_dataset import ParquetDataset

# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.
app = Flask(__name__, static_folder='cirro/client/')
app.register_blueprint(blueprint, url_prefix='/api')


@app.route('/')
def root():
    return send_from_directory(os.path.abspath(os.path.join(app.root_path, "client")), "index.html")


pq_ds = ParquetDataset()
dataset_api.add(pq_ds)
dataset_api.default_provider = pq_ds
auth_api.provider = GoogleAuth()
database_api.provider = FirestoreDatastore()

if __name__ == '__main__':  # for running locally
    from flask_cors import CORS
    CORS(app)
    app.run(host='127.0.0.1', port=5000, debug=True)

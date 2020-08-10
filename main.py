import json
import os

from cirrocumulus.api import blueprint, auth_api, database_api, dataset_api
from cirrocumulus.firestore_datastore import FirestoreDatastore
from cirrocumulus.google_auth import GoogleAuth
from flask import Flask, send_from_directory

# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.


app = Flask(__name__, static_folder='cirrocumulus/client/', static_url_path='')
app.register_blueprint(blueprint, url_prefix='/api')


@app.route('/')
def root():
    return send_from_directory(os.path.abspath(os.path.join(app.root_path, "cirrocumulus", "client")), "index.html")


try:
    from cirrocumulus.parquet_dataset import ParquetDataset

    dataset_api.add(ParquetDataset())
except ModuleNotFoundError:
    pass
try:
    from cirrocumulus.zarr_dataset_backed import ZarrDatasetBacked

    dataset_api.add(ZarrDatasetBacked())
except ModuleNotFoundError:
    pass

with open('cirrocumulus/config.json', 'r') as f:
    config = json.load(f)
auth_api.provider = GoogleAuth(config['clientId'])
database_api.provider = FirestoreDatastore()

if __name__ == '__main__':
    # from flask_cors import CORS
    # CORS(app)
    app.run(host='127.0.0.1', port=5000, debug=True)

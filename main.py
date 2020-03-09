import os

from flask import Flask, send_from_directory

from cirro.api import blueprint, auth_api, database_api, dataset_api
from cirro.firestore_datastore import FirestoreDatastore
from cirro.google_auth import GoogleAuth

# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.


app = Flask(__name__, static_folder='cirro/client/')
app.register_blueprint(blueprint, url_prefix='/api')


@app.route('/')
def root():
    return send_from_directory(os.path.abspath(os.path.join(app.root_path, "cirro", "client")), "index.html")


try:
    from cirro.parquet_dataset import ParquetDataset

    dataset_api.add(ParquetDataset())
except ModuleNotFoundError:
    pass
try:
    from cirro.zarr_dataset import ZarrDataset

    dataset_api.add(ZarrDataset())
except ModuleNotFoundError:
    pass
auth_api.provider = GoogleAuth()
database_api.provider = FirestoreDatastore()

if __name__ == '__main__':
    from flask_cors import CORS

    CORS(app)

    app.run(host='127.0.0.1', port=8080, debug=True)

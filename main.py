from flask import Flask, render_template

from cirro.api import blueprint, auth_api, database_api, dataset_api
from cirro.parquet_backend import ParquetBackend

# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.
app = Flask(__name__, static_folder='build/')
app.register_blueprint(blueprint, url_prefix='/api')


@app.route('/index.html')
@app.route('/')
def root():
    return render_template('index.html')


if __name__ == '__main__':  # for running locally
    app.run(host='127.0.0.1', port=5000, debug=True)
else:
    from cirro.google_auth import GoogleAuth
    from cirro.datastore_api import DatastoreAPI

    dataset_api.add(['pq', 'parquet'], ParquetBackend())
    auth_api.provider = GoogleAuth()
    database_api.provider = DatastoreAPI()

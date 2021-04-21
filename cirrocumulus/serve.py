import os

from cirrocumulus.envir import CIRRO_AUTH_CLIENT_ID, CIRRO_DB_URI, CIRRO_DATABASE, CIRRO_EMAIL, CIRRO_SERVE, \
    CIRRO_FOOTER, CIRRO_UPLOAD, CIRRO_BRAND
from cirrocumulus.launch import create_app

app = None

DEFAULT_DB_URI = 'mongodb://localhost:27017/'
DEFAULT_DATABASE = 'cirrocumulus'


def cached_app():
    global app
    if app is None:
        app = create_app()
        # from flask_cors import CORS
        # CORS(app)
        configure_app(app)
    return app


def configure_app(app):
    from cirrocumulus.api import dataset_api
    from cirrocumulus.no_auth import NoAuth
    auth_client_id = os.environ.get(CIRRO_AUTH_CLIENT_ID)
    db_uri = os.environ.get(CIRRO_DB_URI, DEFAULT_DB_URI)
    database = os.environ.get(CIRRO_DATABASE, DEFAULT_DATABASE)
    email = os.environ.get(CIRRO_EMAIL)
    os.environ[CIRRO_SERVE] = 'true'
    if auth_client_id is None:
        app.config['AUTH'] = NoAuth()
    else:
        from cirrocumulus.google_auth import GoogleAuth
        app.config['AUTH'] = GoogleAuth(auth_client_id)
    from cirrocumulus.mongo_db import MongoDb
    app.config['DATABASE'] = MongoDb(db_uri, database, email)
    from cirrocumulus.parquet_dataset import ParquetDataset
    dataset_api.add(ParquetDataset())
    from cirrocumulus.anndata_dataset import AnndataDataset
    anndata_dataset = AnndataDataset('r' if False else None)
    dataset_api.add(anndata_dataset)


def main(argsv):
    import argparse
    import os
    parser = argparse.ArgumentParser(description='Run cirrocumulus server')
    parser.add_argument('--database', help='MongoDB database', default=DEFAULT_DATABASE)
    parser.add_argument('--db_uri', help='MongoDB database connection URI', default=DEFAULT_DB_URI)
    parser.add_argument('--email', help='Email address that server runs as')
    parser.add_argument('--auth_client_id', help='OAuth client id')
    parser.add_argument('-w', '--workers', dest='workers', help='The number of worker processes', type=int)
    parser.add_argument('-t', '--timeout', dest='timeout',
        help='Workers silent for more than this many seconds are killed and restarted', type=int, default=30)
    parser.add_argument('-b', '--bind', dest='bind',
        help='Server socket to bind. Server sockets can be any of $(HOST), $(HOST):$(PORT), fd://$(FD), or unix:$(PATH). An IP is a valid $(HOST).',
        default='127.0.0.1:5000')
    parser.add_argument('--footer', help='Markdown file to customize the application footer')
    parser.add_argument('--header', help='Markdown file to customize the application header')
    parser.add_argument('--upload', help='URL to allow users to upload files')

    args = parser.parse_args(argsv)

    bind = args.bind if args.bind is not None else '127.0.0.1:5000'
    if args.auth_client_id is not None:
        os.environ[CIRRO_AUTH_CLIENT_ID] = args.auth_client_id

    if args.footer is not None:
        os.environ[CIRRO_FOOTER] = args.footer

    if args.header is not None:
        os.environ[CIRRO_BRAND] = args.header

    if args.email is not None:
        os.environ[CIRRO_EMAIL] = args.email
    if args.db_uri is not None:
        os.environ[CIRRO_DB_URI] = args.db_uri
    if args.database is not None:
        os.environ[CIRRO_DATABASE] = args.database
    if args.workers is not None:
        workers = args.workers
    else:
        import os
        workers = 2 * os.cpu_count()
    if args.upload is not None:
        os.environ[CIRRO_UPLOAD] = args.upload

    run_args = [
            'gunicorn',
            '-b', bind,
            '-w', str(workers),
            '-t', str(args.timeout),
            '-n', 'cirrocumulus-webserver',
            'cirrocumulus.serve:cached_app()'
    ]
    import subprocess
    subprocess.check_call(run_args)


if __name__ == "__main__":
    main()

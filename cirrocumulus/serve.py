
from cirrocumulus.launch import create_app


app = None
db_uri = 'mongodb://localhost:27017/'
database = 'cirrocumulus'
auth_client_id = None
email = None


def cached_app():
    global app
    if not app:
        app = create_app()
        # from flask_cors import CORS
        # CORS(app)
        configure()
    return app


def configure():
    from cirrocumulus.api import dataset_api
    from cirrocumulus.api import auth_api, database_api
    from cirrocumulus.no_auth import NoAuth
    if auth_client_id is None:
        auth_api.provider = NoAuth()
    else:
        from cirrocumulus.google_auth import GoogleAuth
        auth_api.provider = GoogleAuth(auth_client_id)
    from cirrocumulus.mongo_db import MongoDb
    database_api.provider = MongoDb(db_uri, database, email)
    try:
        from cirrocumulus.parquet_dataset import ParquetDataset
        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    try:
        from cirrocumulus.anndata_dataset import AnndataDataset
        anndata_dataset = AnndataDataset('r' if False else None)
        dataset_api.add(anndata_dataset)
    except ModuleNotFoundError:
        pass


def main(argsv):
    global auth_client_id
    global email
    import argparse
    global db_uri
    global database

    parser = argparse.ArgumentParser(description='Run cirrocumulus server')
    parser.add_argument('--database', help='Database')
    parser.add_argument('--db_uri', help='Database connection URI')
    parser.add_argument('--email', help='Email address that server runs as')
    parser.add_argument('--auth_client_id', help='OAuth client id')
    parser.add_argument('-w', '--workers', dest='workers', help='The number of worker processes', type=int)
    parser.add_argument('-b', '--bind', dest='bind',
        help='Server socket to bind. Server sockets can be any of $(HOST), $(HOST):$(PORT), fd://$(FD), or unix:$(PATH). An IP is a valid $(HOST).')

    args = parser.parse_args(argsv)
    bind = args.bind if args.bind is not None else '127.0.0.1:5000'
    if args.auth_client_id is not None:
        auth_client_id = args.auth_client_id
    if args.email is not None:
        email = args.email

    if args.db_uri is not None:
        db_uri = args.db_uri
    if args.database is not None:
        database = args.database

    if args.workers is not None:
        workers = args.workers
    else:
        import os
        workers = 2 * os.cpu_count()
    run_args = [
            'gunicorn',
            '-b', bind,
            '-w', str(workers),
            '-n', 'cirrocumulus-webserver',
            'cirrocumulus.serve:cached_app()'
    ]
    import subprocess
    subprocess.check_call(run_args)


if __name__ == "__main__":
    main()

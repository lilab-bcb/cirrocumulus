from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.launch import create_app
from cirrocumulus.mongo_db import MongoDb

app = None


def cached_app():
    global app
    if not app:
        app = create_app()
        # CORS(app)
        configure()
    return app


def configure():
    from cirrocumulus.api import dataset_api
    from cirrocumulus.api import auth_api, database_api
    from cirrocumulus.no_auth import NoAuth
    auth_api.provider = NoAuth()
    database_api.provider = MongoDb()
    try:
        from cirrocumulus.parquet_dataset import ParquetDataset
        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    try:
        anndataDataset = AnndataDataset('r' if False else None)
        dataset_api.add(anndataDataset)
    except ModuleNotFoundError:
        pass


def main(argsv):
    import argparse
    parser = argparse.ArgumentParser(description='Run cirrocumulus server')
    # parser.add_argument('config', help='Path to config', nargs='+')

    parser.add_argument('--host',
        help='Host IP address')  # set to 0.0.0.0 to make it accessible from other computers WITHOUT SECURITY.

    parser.add_argument('--port', help='Server port', default=5000, type=int)
    # parser.add_argument('--num_workers', help='Number of workers', default=7, type=int)

    args = parser.parse_args(argsv)
    host = args.host if args.host is not None else '127.0.0.1'

    run_args = [
            'gunicorn',
            '-b', host + ':' + str(args.port),
            '-n', 'cirrocumulus-webserver',
            'cirrocumulus.serve:cached_app()'
    ]
    import subprocess
    subprocess.check_call(run_args)


if __name__ == "__main__":
    main()

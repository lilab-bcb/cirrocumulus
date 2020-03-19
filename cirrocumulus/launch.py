from cirrocumulus.anndata_dataset import AnndataDataset


def main(argsv):
    from flask import Flask, send_from_directory
    import argparse
    import webbrowser
    from cirrocumulus.api import blueprint, auth_api, database_api
    from flask_compress import Compress
    app = Flask(__name__, static_folder='client/')
    app.register_blueprint(blueprint, url_prefix='/api')

    @app.route('/')
    def root():
        return send_from_directory(os.path.abspath(os.path.join(app.root_path, "client")), "index.html")

    parser = argparse.ArgumentParser(
        description='Run cirrocumulus')
    parser.add_argument('dataset', help='Path to an h5ad file or parquet file')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--host', help='Host IP address')
    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--no-open', dest='no_open', help='Do not open your web browser', action='store_true')

    args = parser.parse_args(argsv)

    # from flask_cors import CORS
    # CORS(app)
    Compress(app)
    from cirrocumulus.api import dataset_api
    from cirrocumulus.local_db_api import LocalDbAPI
    from cirrocumulus.no_auth import NoAuth
    import os

    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    auth_api.provider = NoAuth()
    database_api.provider = LocalDbAPI(os.path.normpath(args.dataset))
    dataset_api.add(AnndataDataset('r' if args.backed else None))
    try:
        from cirrocumulus.parquet_dataset import ParquetDataset
        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    # dataset_api.add(ZarrDataset())


    if not args.no_open:
        host = args.host if args.host is not None else '127.0.0.1'
        webbrowser.open(host + ':' + str(args.port))
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()

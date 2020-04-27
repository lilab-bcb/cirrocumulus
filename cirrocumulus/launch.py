from cirrocumulus.anndata_dataset import AnndataDataset


def create_app(dataset, backed):
    from cirrocumulus.api import blueprint, auth_api, database_api
    from flask_compress import Compress
    from flask import Flask, send_from_directory
    app = Flask(__name__, static_folder='client/')
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
    app.register_blueprint(blueprint, url_prefix='/api')

    @app.route('/')
    def root():
        return send_from_directory(os.path.abspath(os.path.join(app.root_path, "client")), "index.html")

    # from flask_cors import CORS
    # CORS(app)
    Compress(app)
    from cirrocumulus.api import dataset_api
    from cirrocumulus.local_db_api import LocalDbAPI
    from cirrocumulus.no_auth import NoAuth
    import os

    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    auth_api.provider = NoAuth()
    database_api.provider = LocalDbAPI(os.path.normpath(dataset))

    try:
        from cirrocumulus.parquet_dataset import ParquetDataset
        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    dataset_api.add(AnndataDataset('r' if backed else None))
    return app


def main(argsv):
    import argparse
    import webbrowser

    parser = argparse.ArgumentParser(
        description='Run cirrocumulus')
    parser.add_argument('dataset', help='Path to an h5ad file or parquet file')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--host', help='Host IP address')
    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--processes', help='Number of processes', default=7, type=int)
    parser.add_argument('--no-open', dest='no_open', help='Do not open your web browser', action='store_true')

    args = parser.parse_args(argsv)

    if not args.no_open:
        host = args.host if args.host is not None else '127.0.0.1'
        webbrowser.open(host + ':' + str(args.port))
    app = create_app(args.dataset, args.backed)

    app.run(host=args.host, port=args.port, debug=False, threaded=False, processes=args.processes)
    # from gevent.pywsgi import WSGIServer
    # http_server = WSGIServer(('', args.port), app)
    # http_server.serve_forever()


if __name__ == "__main__":
    main()

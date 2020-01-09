def main(argsv):
    from flask import Flask, send_from_directory
    import argparse
    import webbrowser
    from cirro.api import blueprint, auth_api, database_api
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
    parser.add_argument('--host', help='Host IP address', default="127.0.0.1")
    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--debug', help='Run in debug mode', action='store_true')
    parser.add_argument('--no-open', dest='no_open', help='Do not open your web browser', action='store_true')
    args = parser.parse_args(argsv)

    # from flask_cors import CORS
    # CORS(app)
    Compress(app)
    from cirro.api import dataset_api
    from cirro.h5ad_dataset import H5ADDataset
    from cirro.local_db_api import LocalDbAPI
    from cirro.no_auth import NoAuth
    import os
    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    auth_api.provider = NoAuth()
    database_api.provider = LocalDbAPI([args.dataset])

    try:
        from cirro.parquet_dataset import ParquetDataset
        pq = ParquetDataset()
        dataset_api.add(pq)
        dataset_api.default_provider = pq
    except ModuleNotFoundError:
        pass
    dataset_api.add(H5ADDataset('r' if args.backed else None))
    if not args.no_open:
        url = args.host + ':' + str(args.port)
        webbrowser.open(url)
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()

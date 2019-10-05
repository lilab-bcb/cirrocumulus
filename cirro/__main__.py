def main():
    from flask import Flask, send_from_directory
    import argparse
    from cirro.api import blueprint, auth_api, database_api
    app = Flask(__name__, static_folder='../client/')
    app.register_blueprint(blueprint, url_prefix='/api')

    @app.route('/')
    def root():
        return send_from_directory(os.path.join(app.root_path, "../client"), "index.html")

    parser = argparse.ArgumentParser(
        description='Run cirrocumulus locally')
    parser.add_argument('dataset', help='Path to an h5ad file', action='append')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--host', help='Host IP address', default="127.0.0.1")
    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--debug', help='Run in debug mode', action='store_true')
    args = parser.parse_args()

    # from flask_cors import CORS
    # CORS(app)
    from cirro.api import dataset_api
    from cirro.h5ad_backend import H5ADBackend
    from cirro.local_db_api import LocalDbAPI
    from cirro.no_auth import NoAuth
    import os
    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    auth_api.provider = NoAuth()
    database_api.provider = LocalDbAPI(args.dataset)
    dataset_api.add(['h5ad'], H5ADBackend('r' if args.backed else None))
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()

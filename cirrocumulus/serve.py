import os
import argparse

from cirrocumulus.envir import (
    CIRRO_AUTH,
    CIRRO_AUTH_CLIENT_ID,
    CIRRO_AUTH_PROVIDER,
    CIRRO_BRAND,
    CIRRO_CELL_ONTOLOGY,
    CIRRO_DATABASE,
    CIRRO_DATABASE_CLASS,
    CIRRO_DATASET_PROVIDERS,
    CIRRO_DB_URI,
    CIRRO_FOOTER,
    CIRRO_JOB_RESULTS,
    CIRRO_JOB_TYPE,
    CIRRO_SERVE,
    CIRRO_UPLOAD,
)
from cirrocumulus.launch import create_app
from cirrocumulus.util import add_dataset_providers, create_instance, get_fs


app = None

DEFAULT_DB_URI = "mongodb://localhost:27017/cirrocumulus"


def cached_app():
    global app
    if app is None:
        app = create_app()
        configure_app(app)
    return app


def configure_app(app):
    from cirrocumulus.no_auth import NoAuth

    auth_client_id = os.environ.get(CIRRO_AUTH_CLIENT_ID)

    os.environ[CIRRO_SERVE] = "true"
    os.environ[CIRRO_JOB_TYPE + "de"] = "cirrocumulus.job_api.run_de"
    if auth_client_id is None:
        app.config[CIRRO_AUTH] = NoAuth()
    else:
        if os.environ.get(CIRRO_AUTH_PROVIDER, "google").lower() == "google":
            from cirrocumulus.google_auth import GoogleAuth

            app.config[CIRRO_AUTH] = GoogleAuth()
        elif os.environ.get(CIRRO_AUTH_PROVIDER, "").lower() == "okta":
            from cirrocumulus.okta_auth import OktaAuth

            app.config[CIRRO_AUTH] = OktaAuth()
        else:
            raise ValueError(
                "Unknown CIRRO_AUTH_PROVIDER - {}".format(os.environ[CIRRO_AUTH_PROVIDER])
            )
    if os.environ.get(CIRRO_DATABASE_CLASS) is None:
        os.environ[CIRRO_DATABASE_CLASS] = "cirrocumulus.mongo_db.MongoDb"
    if os.environ[CIRRO_DB_URI] == "":
        os.environ[CIRRO_DATABASE_CLASS] = "cirrocumulus.local_db_api.LocalDbAPI"

    app.config[CIRRO_DATABASE] = create_instance(os.environ[CIRRO_DATABASE_CLASS])
    os.environ[CIRRO_DATASET_PROVIDERS] = ",".join(
        [
            "cirrocumulus.parquet_dataset.ParquetDataset",
            "cirrocumulus.zarr_dataset.ZarrDataset",
            "cirrocumulus.tiledb_dataset.TileDBDataset",
        ]
    )
    add_dataset_providers()


def create_parser(description=False):
    parser = argparse.ArgumentParser(description="Run cirrocumulus server" if description else None)
    parser.add_argument("--db_uri", help="Database connection URI", default=DEFAULT_DB_URI)
    parser.add_argument(
        "-w", "--workers", dest="workers", help="The number of worker processes", type=int
    )
    parser.add_argument(
        "-t",
        "--timeout",
        dest="timeout",
        help="Workers silent for more than this many seconds are killed and restarted",
        type=int,
        default=30,
    )
    parser.add_argument(
        "-b",
        "--bind",
        dest="bind",
        help="Server socket to bind. Server sockets can be any of $(HOST), $(HOST):$(PORT), fd://$(FD), or unix:$(PATH). An IP is a valid $(HOST).",
        default="127.0.0.1:5000",
    )
    parser.add_argument("--footer", help="Markdown file to customize the application footer")
    parser.add_argument("--header", help="Markdown file to customize the application header")
    parser.add_argument("--upload", help="URL to allow users to upload files")
    parser.add_argument(
        "--results", help="URL to save user computed results (e.g. differential expression) to"
    )
    parser.add_argument("--ontology", help="Path to ontology in OBO format for annotation")
    return parser


def main(argsv):
    args = create_parser(True).parse_args(argsv)

    bind = args.bind if args.bind is not None else "127.0.0.1:5000"
    if args.ontology is not None:
        os.environ[CIRRO_CELL_ONTOLOGY] = args.ontology
    os.environ[CIRRO_DB_URI] = args.db_uri

    if args.footer is not None:
        os.environ[CIRRO_FOOTER] = args.footer

    if args.header is not None:
        os.environ[CIRRO_BRAND] = args.header

    if args.workers is not None:
        workers = args.workers
    else:
        workers = 2 * os.cpu_count()
    if args.upload is not None:
        os.environ[CIRRO_UPLOAD] = args.upload
    if args.results is not None:
        os.environ[CIRRO_JOB_RESULTS] = args.results
        get_fs(os.environ[CIRRO_JOB_RESULTS]).makedirs(os.environ[CIRRO_JOB_RESULTS], exist_ok=True)

    run_args = [
        "gunicorn",
        "-b",
        bind,
        "-w",
        str(workers),
        "-t",
        str(args.timeout),
        "-n",
        "cirrocumulus-webserver",
        "cirrocumulus.serve:cached_app()",
    ]
    # if args.gunicorn is not None:
    #     run_args += args.gunicorn.split(' ')
    import subprocess

    subprocess.check_call(run_args)


if __name__ == "__main__":
    main()

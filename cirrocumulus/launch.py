import os
import argparse

import anndata

import cirrocumulus
from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.envir import (
    CIRRO_AUTH,
    CIRRO_CELL_ONTOLOGY,
    CIRRO_COMPRESS,
    CIRRO_DATABASE,
    CIRRO_JOB_RESULTS,
    CIRRO_JOB_TYPE,
)
from cirrocumulus.io_util import SPATIAL_HELP, add_spatial, filter_markers, get_markers
from cirrocumulus.local_db_api import LocalDbAPI
from cirrocumulus.util import get_fs


def configure_app(app, list_of_dataset_paths, spatial_directories, marker_paths):
    from cirrocumulus.api import dataset_api
    from cirrocumulus.no_auth import NoAuth

    try:
        from cirrocumulus.parquet_dataset import ParquetDataset

        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    try:
        from cirrocumulus.zarr_dataset import ZarrDataset

        dataset_api.add(ZarrDataset())
    except ModuleNotFoundError:
        pass
    app.config[CIRRO_AUTH] = NoAuth()
    os.environ[CIRRO_JOB_TYPE + "de"] = "cirrocumulus.job_api.run_de"
    os.environ[CIRRO_JOB_TYPE + "ot_trajectory"] = "cirrocumulus.job_api.run_ot_trajectory"
    anndata_dataset = AnndataDataset()
    dataset_api.add(anndata_dataset)
    dataset_ids = []
    for dataset_paths in list_of_dataset_paths:
        dataset_paths = dataset_paths.split(",")
        dataset_id = dataset_paths[0]
        dataset_ids.append(dataset_id)
        if len(dataset_paths) > 1:
            datasets = []
            for i in range(len(dataset_paths)):
                dataset = anndata_dataset.get_data(get_fs(dataset_paths[i]), dataset_paths[i])
                if "group" not in dataset.var:
                    dataset.var["group"] = dataset.uns.get("name", "dataset {}".format(i + 1))
                datasets.append(dataset)
            adata = anndata.concat(datasets, axis=1, label="group", merge="unique")
            dataset.obsm = datasets[0].obsm
            adata.var.index = adata.var.index.str.replace("/", "_")
            adata.var_names_make_unique()
            anndata_dataset.add_data(dataset_id, adata)

    app.config[CIRRO_DATABASE] = LocalDbAPI(dataset_ids)

    if spatial_directories is not None and len(spatial_directories) > 0:
        for i in range(len(spatial_directories)):
            spatial_directory = spatial_directories[i]
            if spatial_directory != "":
                adata = anndata_dataset.get_data(get_fs(dataset_ids[i]), dataset_ids[i])
                if not add_spatial(adata, spatial_directory):
                    print("No spatial data found in {}".format(spatial_directory))

    if marker_paths is not None and len(marker_paths) > 0:
        markers = get_markers(marker_paths)
        for dataset_id in dataset_ids:
            d = anndata_dataset.get_data(get_fs(dataset_id), dataset_id)
            existing_markers = d.uns.get("markers", [])
            markers += existing_markers
            # remove genes in dict that are not in dataset
            d.uns["markers"] = filter_markers(d, markers)


def create_app():
    from flask import Flask, send_from_directory
    from flask_compress import Compress

    from cirrocumulus.api import cirro_blueprint

    app = Flask(
        __name__, static_folder=os.path.join(cirrocumulus.__path__[0], "client"), static_url_path=""
    )
    app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0
    app.register_blueprint(cirro_blueprint, url_prefix="/api")

    @app.route("/")
    def root():
        return send_from_directory(
            os.path.abspath(os.path.join(app.root_path, "client")), "index.html"
        )

    if bool(os.environ.get(CIRRO_COMPRESS, "true")):
        Compress(app)
    return app


def create_parser(description=False):
    parser = argparse.ArgumentParser(description="Run cirrocumulus" if description else None)
    parser.add_argument(
        "dataset",
        help="Path(s) to dataset in h5ad, loom, Seurat, TileDB, zarr, or STAR-Fusion format. Separate multiple datasets with "
        "a comma instead of a space in order to join datasets by cell id",
        nargs="+",
    )
    parser.add_argument("--spatial", help=SPATIAL_HELP, nargs="*")
    parser.add_argument(
        "--markers",
        help='Path(s) to JSON file that maps name to features. For example {"a":["gene1", "gene2"], "b":["gene3"]}',
        nargs="*",
    )
    parser.add_argument(
        "--host", help="Host IP address"
    )  # set to 0.0.0.0 to make it accessible from other computers WITHOUT SECURITY.

    parser.add_argument("--port", help="Server port", default=5000, type=int)
    parser.add_argument(
        "--no-open", dest="no_open", help="Do not open your web browser", action="store_true"
    )
    parser.add_argument(
        "--results", help="URL to save user computed results (e.g. differential expression)"
    )
    parser.add_argument("--ontology", help="Path to ontology in OBO format for annotation")
    parser.add_argument(
        "--tmap", help="Path(s) to transport maps directory computed with WOT", nargs="*"
    )
    return parser


def main(argsv):
    args = create_parser(True).parse_args(argsv)
    if args.results is not None:
        os.environ[CIRRO_JOB_RESULTS] = args.results
    else:
        os.environ[CIRRO_JOB_RESULTS] = os.path.join(
            os.path.dirname(args.dataset[0].rstrip("/")), "results"
        )
    get_fs(os.environ[CIRRO_JOB_RESULTS]).makedirs(os.environ[CIRRO_JOB_RESULTS], exist_ok=True)
    if args.ontology is not None:
        os.environ[CIRRO_CELL_ONTOLOGY] = args.ontology
    app = create_app()
    configure_app(app, args.dataset, args.spatial, args.markers)
    if args.tmap is not None:
        tmaps = []
        unique_names = set()
        for tmap_dir in args.tmap:
            name = os.path.basename(tmap_dir.rstrip("/"))
            counter = 1
            while name in unique_names:
                name = f"{name}-{counter}"
                counter = counter + 1
            tmaps.append({"path": tmap_dir, "name": name})
        app.config[CIRRO_DATABASE].add_metadata(None, {"ot": {"tmaps": tmaps}})
    if not args.no_open:
        import webbrowser

        import requests

        host = args.host if args.host is not None else "http://127.0.0.1"
        url = host + ":" + str(args.port)
        try:
            if requests.get(url).ok:
                import sys

                sys.exit("Address already in use")
        except:
            pass
        webbrowser.open(url)
    from flask import cli

    cli.show_server_banner = lambda *_: None  # suppress warning message
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])

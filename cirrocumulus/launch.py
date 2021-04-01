import os

from cirrocumulus.anndata_dataset import AnndataDataset
from cirrocumulus.io_util import get_markers, filter_markers, add_spatial, SPATIAL_HELP
from cirrocumulus.parquet_dataset import ParquetDataset


def configure_app(app, list_of_dataset_paths, spatial_directories, backed, marker_paths):
    from cirrocumulus.api import dataset_api
    from cirrocumulus.local_db_api import LocalDbAPI
    from cirrocumulus.no_auth import NoAuth

    app.config['AUTH'] = NoAuth()
    dataset_api.add(ParquetDataset())
    anndata_dataset = AnndataDataset('r' if backed else None)
    dataset_ids = []
    for dataset_paths in list_of_dataset_paths:
        dataset_paths = dataset_paths.split(',')
        dataset_id = os.path.normpath(dataset_paths[0])
        dataset_ids.append(dataset_id)
        if len(dataset_paths) > 1:
            to_concat = []
            all_ids = None
            for path in dataset_paths:
                d = anndata_dataset.get_data(path)
                all_ids = d.obs.index.union(all_ids) if all_ids is not None else d.obs.index
                to_concat.append(d)
            for i in range(len(to_concat)):
                d = to_concat[i]
                missing_ids = all_ids.difference(d.obs.index)
                if len(missing_ids) > 0:
                    import scipy.sparse
                    import anndata
                    import pandas as pd
                    X = None
                    if d.shape[1] > 0:
                        empty = scipy.sparse.csr_matrix((len(missing_ids), d.shape[1]))
                        X = scipy.sparse.vstack((d.X, empty), format='csr')
                    missing_df = pd.DataFrame(index=missing_ids)
                    for column in d.obs:
                        if pd.api.types.is_bool_dtype(d.obs[column]):
                            missing_df[column] = False

                    obs = pd.concat((d.obs, missing_df))
                    # for column in d.obs:
                    #     if pd.api.types.is_categorical_dtype(d.obs[column]):
                    #         obs[column] = obs[column].astype('category')
                    d = anndata.AnnData(X=X, obs=obs, var=d.var)
                d = d[all_ids]  # same order
                to_concat[i] = d
            X_list = []
            obs = None
            obsm = {}
            var = None
            for d in to_concat:
                if d.shape[1] > 0:
                    X_list.append(d.X)
                    var = pd.concat((var, d.var)) if var is not None else d.var
                obs = obs.join(d.obs) if obs is not None else d.obs
                for key in d.obsm_keys():
                    obsm[key] = d.obsm[key]

            X = scipy.sparse.hstack(X, format='csr') if len(X_list) > 1 else X_list[0]
            adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm)
            adata.var_names_make_unique()
            anndata_dataset.add_data(dataset_id, adata)
        dataset_api.add(anndata_dataset)

    app.config['DATABASE'] = LocalDbAPI(dataset_ids)

    if spatial_directories is not None and len(spatial_directories) > 0:
        for i in range(len(spatial_directories)):
            spatial_directory = spatial_directories[i]
            if spatial_directory != '':
                adata = anndata_dataset.get_data(dataset_ids[i])
                if not add_spatial(adata, spatial_directory):
                    print('No spatial data found in {}'.format(spatial_directory))

    if marker_paths is not None and len(marker_paths) > 0:
        markers = get_markers(marker_paths)
        for dataset_id in dataset_ids:
            d = anndata_dataset.get_data(dataset_id)
            existing_markers = d.uns.get('markers', [])
            markers += existing_markers
            # remove genes in dict that are not in dataset
            d.uns['markers'] = filter_markers(d, markers)


def create_app():
    from cirrocumulus.api import blueprint
    from flask_compress import Compress
    from flask import Flask, send_from_directory
    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    app = Flask(__name__, static_folder='client', static_url_path='')
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
    app.register_blueprint(blueprint, url_prefix='/api')

    @app.route('/')
    def root():
        return send_from_directory(os.path.abspath(os.path.join(app.root_path, "client")), "index.html")

    Compress(app)
    return app


def main(argsv):
    import argparse
    parser = argparse.ArgumentParser(description='Run cirrocumulus')
    parser.add_argument('dataset',
        help='Path to dataset in h5ad, loom, or STAR-Fusion format. Separate multiple datasets with '
             'a comma instead of a space in order to join datasets by cell id', nargs='+')
    parser.add_argument('--spatial', help=SPATIAL_HELP, nargs='*')
    parser.add_argument('--markers',
        help='Path to JSON file that maps name to features. For example {"a":["gene1", "gene2"], "b":["gene3"]}',
        nargs='*')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--host',
        help='Host IP address')  # set to 0.0.0.0 to make it accessible from other computers WITHOUT SECURITY.

    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--no-open', dest='no_open', help='Do not open your web browser', action='store_true')

    args = parser.parse_args(argsv)
    app = create_app()
    # from flask_cors import CORS
    # CORS(app)
    configure_app(app, args.dataset, args.spatial, args.backed, args.markers)
    if not args.no_open:
        import webbrowser
        host = args.host if args.host is not None else 'http://127.0.0.1'
        url = host + ':' + str(args.port)
        webbrowser.open(url)
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()

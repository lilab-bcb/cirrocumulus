from cirrocumulus.anndata_dataset import AnndataDataset


def create_app(dataset_paths, backed):
    from cirrocumulus.api import blueprint, auth_api, database_api
    from flask_compress import Compress
    from flask import Flask, send_from_directory
    from cirrocumulus.api import dataset_api
    from cirrocumulus.local_db_api import LocalDbAPI
    from cirrocumulus.no_auth import NoAuth
    import os

    os.environ['WERKZEUG_RUN_MAIN'] = 'true'
    auth_api.provider = NoAuth()
    database_api.provider = LocalDbAPI(os.path.normpath(dataset_paths[0]))
    try:
        from cirrocumulus.parquet_dataset import ParquetDataset
        dataset_api.add(ParquetDataset())
    except ModuleNotFoundError:
        pass
    anndataDataset = AnndataDataset('r' if backed else None)
    if len(dataset_paths) > 1:
        to_concat = []
        all_ids = None
        for path in dataset_paths:
            d = anndataDataset.get_data(path)
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
        anndataDataset.add_data(os.path.normpath(dataset_paths[0]), adata)
    dataset_api.add(anndataDataset)

    app = Flask(__name__, static_folder='client', static_url_path='')
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
    app.register_blueprint(blueprint, url_prefix='/api')

    @app.route('/')
    def root():
        return send_from_directory(os.path.abspath(os.path.join(app.root_path, "client")), "index.html")

    # from flask_cors import CORS
    # CORS(app)
    Compress(app)
    return app


def main(argsv):
    import argparse
    parser = argparse.ArgumentParser(description='Run cirrocumulus')
    parser.add_argument('dataset', help='Path to dataset', nargs='+')
    parser.add_argument('--backed', help='Load h5ad file in backed mode', action='store_true')
    parser.add_argument('--host', help='Host IP address')
    parser.add_argument('--port', help='Server port', default=5000, type=int)
    parser.add_argument('--processes', help='Number of processes', default=7, type=int)
    parser.add_argument('--no-open', dest='no_open', help='Do not open your web browser', action='store_true')

    args = parser.parse_args(argsv)
    app = create_app(args.dataset, args.backed)
    if not args.no_open:
        import webbrowser
        host = args.host if args.host is not None else 'http://127.0.0.1'
        url = host + ':' + str(args.port)
        webbrowser.open(url)
    app.run(host=args.host, port=args.port, debug=False, threaded=False, processes=args.processes)


if __name__ == "__main__":
    main()

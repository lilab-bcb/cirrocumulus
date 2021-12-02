import argparse
import logging
import os

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
from cirrocumulus.anndata_util import get_scanpy_marker_keys, dataset_schema, ADATA_MODULE_UNS_KEY
from cirrocumulus.io_util import get_markers, filter_markers, add_spatial, SPATIAL_HELP, unique_id
from cirrocumulus.util import to_json, get_fs, open_file

logger = logging.getLogger("cirro")

cluster_fields = ['anno', 'cell_type', 'celltype', 'leiden', 'louvain', 'seurat_cluster', 'cluster']
categorical_fields_convert = ['seurat_clusters']


def read_adata(path, spatial_directory=None, use_raw=False):
    if path.lower().endswith('.loom'):
        adata = anndata.read_loom(path)
    elif path.lower().endswith('.zarr'):
        adata = anndata.read_zarr(path)
    else:
        adata = anndata.read(path)
    if 'module' in adata.uns:
        adata.uns[ADATA_MODULE_UNS_KEY] = anndata.AnnData(X=adata.uns['module']['X'],
                                                          var=adata.uns['module']['var'])
    if use_raw and adata.raw is not None and adata.shape[0] == adata.raw.shape[0]:
        logger.info('Using adata.raw')
        adata = anndata.AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs, obsm=adata.obsm, uns=adata.uns)

    if spatial_directory is not None:
        if not add_spatial(adata, spatial_directory):
            logger.info('No spatial data found in {}'.format(spatial_directory))

    for field in categorical_fields_convert:
        if field in adata.obs and not pd.api.types.is_categorical_dtype(adata.obs[field]):
            logger.info('Converting {} to categorical'.format(field))
            adata.obs[field] = adata.obs[field].astype(str).astype('category')
    return adata


class PrepareData:

    def __init__(self, datasets, output, dimensions=None, groups=[], group_nfeatures=10, markers=[],
                 output_format='parquet', no_auto_groups=False, save_whitelist=None):
        self.groups = groups
        self.group_nfeatures = group_nfeatures
        self.markers = markers
        self.output_format = output_format
        self.no_auto_groups = no_auto_groups
        self.save_whitelist = save_whitelist

        for dataset in datasets:
            for key in list(dataset.obsm.keys()):
                m = dataset.obsm[key]
                dim = m.shape[1]
                if not (1 < dim <= 3):
                    del dataset.obsm[key]
        if len(datasets) > 1:
            for i in range(len(datasets)):
                dataset = datasets[i]
                if 'group' not in dataset.var:
                    dataset.var['group'] = dataset.uns.get('name', 'dataset {}'.format(i + 1))
                if i > 0 and not np.array_equal(datasets[0].obs.index, dataset.obs.index):
                    raise ValueError('obs ids are not equal')
            dataset = anndata.concat(datasets, axis=1, label='group', merge='unique')
        else:
            dataset = datasets[0]
        dataset.var.index = dataset.var.index.str.replace('/', '_')
        dataset.var_names_make_unique()
        dataset.obs.index.name = 'id'
        dataset.var.index.name = 'id'
        self.base_output = output
        dimensions_supplied = dimensions is not None and len(dimensions) > 0
        self.dimensions = [] if not dimensions_supplied else dimensions
        self.measures = []
        self.others = []
        self.dataset = dataset
        if scipy.sparse.issparse(dataset.X) and not scipy.sparse.isspmatrix_csc(dataset.X):
            dataset.X = dataset.X.tocsc()
        for i in range(len(dataset.obs.columns)):
            name = dataset.obs.columns[i]
            c = dataset.obs[name]
            if pd.api.types.is_object_dtype(c):
                dataset.obs[name] = dataset.obs[name].astype('category')
                c = dataset.obs[name]
            if not dimensions_supplied and pd.api.types.is_categorical_dtype(c):
                if 1 < len(c.cat.categories) < 2000:
                    self.dimensions.append(name)
                    if c.isna().sum() > 0:
                        logger.info('Replacing nans in {}'.format(name))
                        dataset.obs[name] = dataset.obs[name].astype(str)
                        dataset.obs.loc[dataset.obs[name].isna(), name] = ''
                        dataset.obs[name] = dataset.obs[name].astype('category')
                else:
                    self.others.append(name)
            elif not pd.api.types.is_string_dtype(c) and not pd.api.types.is_object_dtype(c):
                self.measures.append('obs/' + name)
            else:
                self.others.append(name)

    def execute(self):
        output_format = self.output_format
        dataset = self.dataset
        if self.groups is None and not self.no_auto_groups:
            groups = []
            existing_fields = set()
            scanpy_marker_keys = get_scanpy_marker_keys(dataset)
            for key in scanpy_marker_keys:
                group_by = dataset.uns[key]['params']['groupby']
                if isinstance(group_by, np.ndarray):
                    group_by = ','.join(group_by)
                existing_fields.add(group_by)
            for field in dataset.obs.columns:
                field_lc = field.lower()
                for cluster_field in cluster_fields:
                    if field_lc.find(cluster_field) != -1 and cluster_field not in existing_fields:
                        groups.append(field)
                        break
            self.groups = groups
        if self.groups is not None and len(self.groups) > 0:
            use_pegasus = False
            use_scanpy = False
            try:
                import pegasus as pg
                use_pegasus = True
            except ModuleNotFoundError:
                pass
            if not use_pegasus:
                try:
                    import scanpy as sc
                    use_scanpy = True
                except ModuleNotFoundError:
                    pass
            if not use_pegasus and not use_scanpy:
                raise ValueError('Please install pegasuspy or scanpy to compute markers')

            for group in self.groups:
                field = group
                if group not in dataset.obs:  # test if multiple comma separated fields
                    split_groups = group.split(',')
                    if len(split_groups) > 1:
                        use_split_groups = True
                        for split_group in split_groups:
                            if split_group not in dataset.obs:
                                use_split_groups = False
                                break
                        if use_split_groups:
                            dataset.obs[field] = dataset.obs[split_groups[0]].str.cat(dataset.obs[split_groups[1:]],
                                                                                      sep=',')

                if field in dataset.obs:
                    if not pd.api.types.is_categorical_dtype(dataset.obs[field]):
                        dataset.obs[field] = dataset.obs[field].astype(str).astype('category')
                    if len(dataset.obs[field].cat.categories) > 1:
                        logger.info('Computing markers for {}'.format(field))
                        key_added = 'rank_genes_' + str(field)
                        if use_pegasus:
                            pg.de_analysis(dataset, cluster=field, de_key=key_added)
                        else:
                            sc.tl.rank_genes_groups(dataset, field, key_added=key_added, method='t-test')
                else:
                    logger.info(group + ' not found in ' + ', '.join(dataset.obs.columns))
        schema = self.get_schema()
        schema['format'] = output_format
        if output_format in ['parquet', 'zarr']:
            output_dir = self.base_output
        else:
            output_dir = os.path.splitext(self.base_output)[0]
        filesystem = get_fs(output_dir)
        filesystem.makedirs(output_dir, exist_ok=True)
        results = schema.get('results', [])

        if len(results) > 0:
            uns_dir = os.path.join(output_dir, 'uns')
            is_gzip = output_format != 'jsonl'
            filesystem.makedirs(uns_dir, exist_ok=True)

            for i in range(len(results)):
                full_result = results[i]
                result_id = full_result.pop('id')
                # keep id, name, type in schema, store rest externally
                results[i] = dict(id=result_id, name=full_result.pop('name'), type=full_result.pop('type'),
                                  content_type='application/json', content_encoding='gzip' if is_gzip else None)
                json_result = to_json(full_result)

                result_path = os.path.join(uns_dir, result_id + '.json.gz') if is_gzip else os.path.join(uns_dir,
                                                                                                             result_id + '.json')
                with open_file(result_path, 'wt', compression='gzip' if is_gzip else None) as out:
                    out.write(json_result)
        images = dataset.uns.get('images')
        if images is not None:
            image_dir = os.path.join(output_dir, 'images')
            filesystem.makedirs(image_dir, exist_ok=True)
            for image in images:
                src = image['image']
                dest = os.path.join(image_dir, os.path.basename(src))
                filesystem.copy(src, dest)
                image['image'] = 'images/' + os.path.basename(src)

        if output_format == 'parquet':
            from cirrocumulus.parquet_output import save_dataset_pq
            save_dataset_pq(dataset, schema, self.base_output, filesystem, self.save_whitelist)
        elif output_format == 'jsonl':
            from cirrocumulus.jsonl_io import save_dataset_jsonl
            save_dataset_jsonl(dataset, schema, output_dir, self.base_output, filesystem)
        elif output_format == 'zarr':
            from cirrocumulus.zarr_output import save_dataset_zarr
            save_dataset_zarr(dataset, schema, self.base_output, filesystem, self.save_whitelist)
        else:
            raise ValueError("Unknown format")

    def get_schema(self):
        result = dataset_schema(self.dataset, n_features=self.group_nfeatures)
        markers = result.get('markers', [])

        if self.markers is not None:  # add results specified from file
            markers += get_markers(self.markers)
            markers = filter_markers(self.dataset, markers)

        for marker in markers:
            if marker.get('id') is None:
                marker['id'] = unique_id()
            marker['readonly'] = True
        result['markers'] = markers
        result['format'] = self.output_format
        return result


def main(argsv):
    parser = argparse.ArgumentParser(
        description='Prepare a dataset for cirrocumulus server.')
    parser.add_argument('dataset', help='Path to a h5ad, loom, or Seurat file', nargs='+')
    parser.add_argument('--out', help='Path to output directory')
    parser.add_argument('--format', help='Output format', choices=['parquet', 'jsonl', 'zarr'],
                        default='zarr')
    parser.add_argument('--whitelist',
                        help='Optional whitelist of fields to save. Only applies when output format is parquet',
                        choices=['obs', 'obsm', 'X'],
                        action='append')
    parser.add_argument('--markers',
                        help='Path to JSON file of precomputed markers that maps name to features. For example {"a":["gene1", "gene2"], "b":["gene3"]',
                        action='append')
    parser.add_argument('--no-auto-groups', dest='no_auto_groups',
                        help='Disable automatic cluster field detection to compute differential expression results for',
                        action='store_true')
    parser.add_argument('--groups',
                        help='List of groups to compute markers for (e.g. louvain). Markers created with cumulus/scanpy are automatically included. Separate multiple groups with a comma to combine groups using "AND" logic (e.g. louvain,day)',
                        action='append')
    parser.add_argument('--group_nfeatures', help='Number of marker genes/features to include', type=int, default=10)
    parser.add_argument('--spatial', help=SPATIAL_HELP)
    args = parser.parse_args(argsv)
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())
    out = args.out
    no_auto_groups = args.no_auto_groups
    save_whitelist = args.whitelist
    input_datasets = args.dataset  # multimodal
    output_format = args.format
    if out is None:
        out = os.path.splitext(os.path.basename(input_datasets[0]))[0]
    if out.endswith('/'):
        out = out[:len(out) - 1]
    output_format2extension = dict(parquet='.cpq', jsonl='.jsonl', zarr='.zarr', h5ad='.h5ad')
    if not out.lower().endswith(output_format2extension[output_format]):
        out += output_format2extension[output_format]

    datasets = []
    tmp_files = []
    for input_dataset in input_datasets:
        use_raw = False
        if input_dataset.lower().endswith('.rds'):
            import subprocess
            import tempfile
            import pkg_resources

            _, h5_file = tempfile.mkstemp(suffix='.h5ad')
            os.remove(h5_file)
            subprocess.check_call(
                ['Rscript', pkg_resources.resource_filename("cirrocumulus", 'seurat2h5ad.R'), input_dataset, h5_file])
            input_dataset = h5_file
            tmp_file = h5_file
            use_raw = True
            tmp_files.append(tmp_file)
        adata = read_adata(input_dataset, spatial_directory=args.spatial, use_raw=use_raw)
        datasets.append(adata)
        adata.uns['name'] = os.path.splitext(os.path.basename(input_dataset.rstrip('/')))[0]

    prepare_data = PrepareData(datasets=datasets, output=out, dimensions=args.groups, groups=args.groups,
                               group_nfeatures=args.group_nfeatures,
                               markers=args.markers, output_format=output_format, no_auto_groups=no_auto_groups,
                               save_whitelist=save_whitelist)
    prepare_data.execute()
    for tmp_file in tmp_files:
        os.remove(tmp_file)


if __name__ == '__main__':
    import sys

    main(sys.argv)

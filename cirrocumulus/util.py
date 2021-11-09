# import pandas.io.json as json
import os
from urllib.parse import urlparse

import fsspec
import numpy as np
import pandas as pd
import pandas._libs.json as ujson
from flask import make_response

from cirrocumulus.envir import CIRRO_DATASET_PROVIDERS


def add_dataset_providers():
    from cirrocumulus.api import dataset_api
    dataset_providers = []

    for p in os.environ[CIRRO_DATASET_PROVIDERS].split(','):
        try:
            dataset_api.add(create_instance(p))
            dataset_providers.append(p)
        except ModuleNotFoundError:  # ignore if required libraries are not installed
            pass

    os.environ[CIRRO_DATASET_PROVIDERS] = ','.join(dataset_providers)


def import_path(name):
    import importlib
    dot_index = name.rfind('.')
    return getattr(importlib.import_module(name[0:dot_index]), name[dot_index + 1:])


def create_instance(class_name):
    return import_path(class_name)()


def get_scheme(path):
    pr = urlparse(path)
    if len(pr.scheme) <= 1:  # for file paths: /foo/bar/test.h5ad or C:/foo/bar/test.h5ad
        return 'file'
    return pr.scheme


fsspec_kwargs = dict(s3_additional_kwargs=dict(ACL='bucket-owner-full-control'))


def get_fs(path):
    return fsspec.filesystem(get_scheme(path), **fsspec_kwargs)


def open_file(urlpath, mode="rb", compression=None):
    return fsspec.open(urlpath, mode=mode, compression=compression, **fsspec_kwargs)


def to_json(data, orient='values'):
    return ujson.dumps(data, double_precision=2, orient=orient)


def json_response(data, response=200):
    # response = make_response(simplejson.dumps(data, check_circular=True), response)
    # response = make_response(json.dumps(data), response)
    s = ujson.dumps(data, double_precision=2, orient='values')
    # s = nujson.dumps(data, double_precision=1)
    response = make_response(s, response)
    response.headers['Content-Type'] = 'application/json'
    return response


def get_email_domain(email):
    at_index = email.find('@')
    domain = None
    if at_index != -1:
        domain = email[at_index + 1:]
    return domain


def load_dataset_schema(url):
    import os

    import json

    def get_extension(path):
        name, ext = os.path.splitext(path)
        if ext == '.gz':
            name, ext = os.path.splitext(name)
            if ext == '.json':
                ext = '.json.gz'
        return ext

    extension = get_extension(url)
    json_schema = None
    if extension in ['.json', '.json.gz', '']:
        scheme = get_scheme(url)
        fs = get_fs(scheme)
        if extension == '':
            url = os.path.join(url, 'index.json.gz')
            extension = get_extension(url)
        if extension == '.json.gz':
            import gzip
            with gzip.open(fs.open(url)) as f:
                json_schema = json.load(f)
        else:
            with fs.open(url) as f:
                json_schema = json.load(f)
    return json_schema


def write_top_half_gct(f, row_metadata_df, col_metadata_df, metadata_null, filler_null):
    """ Write the top half of the gct file: top-left filler values, row metadata
    headers, and top-right column metadata.
    Args:
        f (file handle): handle for output file
        row_metadata_df (pandas df)
        col_metadata_df (pandas df)
        metadata_null (string): how to represent missing values in the metadata
        filler_null (string): what value to fill the top-left filler block with
    Returns:
        None
    """
    # Initialize the top half of the gct including the third line
    size_of_top_half_df = (1 + col_metadata_df.shape[1],
                           1 + row_metadata_df.shape[1] + col_metadata_df.shape[0])

    top_half_df = pd.DataFrame(np.full(size_of_top_half_df, filler_null, dtype=object))

    # Assemble the third line of the gct: "id", then rhds, then cids
    top_half_df.iloc[0, :] = np.hstack(("id", row_metadata_df.columns.values, col_metadata_df.index.values))

    # Insert the chds
    top_half_df.iloc[range(1, top_half_df.shape[0]), 0] = col_metadata_df.columns.values

    # Insert the column metadata, but first convert to strings and replace NaNs
    col_metadata_indices = (range(1, top_half_df.shape[0]),
                            range(1 + row_metadata_df.shape[1], top_half_df.shape[1]))
    # pd.DataFrame.at to insert into dataframe(python3)
    top_half_df.at[col_metadata_indices[0], col_metadata_indices[1]] = (
        col_metadata_df.astype(str).replace("nan", value=metadata_null).T.values)

    # Write top_half_df to file
    top_half_df.to_csv(f, header=False, index=False, sep="\t")


def write_bottom_half_gct(f, row_metadata_df, data_df, data_null, data_float_format, metadata_null):
    """ Write the bottom half of the gct file: row metadata and data.
    Args:
        f (file handle): handle for output file
        row_metadata_df (pandas df)
        data_df (pandas df)
        data_null (string): how to represent missing values in the data
        metadata_null (string): how to represent missing values in the metadata
        data_float_format (string): how many decimal points to keep in representing data
    Returns:
        None
    """
    # create the left side of the bottom half of the gct (for the row metadata)
    size_of_left_bottom_half_df = (row_metadata_df.shape[0],
                                   1 + row_metadata_df.shape[1])
    left_bottom_half_df = pd.DataFrame(np.full(size_of_left_bottom_half_df, metadata_null, dtype=object))

    # create the full bottom half by combining with the above with the matrix data
    bottom_half_df = pd.concat([left_bottom_half_df, data_df.reset_index(drop=True)], axis=1)
    bottom_half_df.columns = range(bottom_half_df.shape[1])

    # Insert the rids
    bottom_half_df.iloc[:, 0] = row_metadata_df.index.values

    # Insert the row metadata, but first convert to strings and replace NaNs
    row_metadata_col_indices = range(1, 1 + row_metadata_df.shape[1])
    bottom_half_df.iloc[:, row_metadata_col_indices] = (
        row_metadata_df.astype(str).replace("nan", value=metadata_null).values)

    # Write bottom_half_df to file
    bottom_half_df.to_csv(f, header=False, index=False, sep="\t",
                          na_rep=data_null,
                          float_format=data_float_format)


def adata2gct(adata, f):
    data_float_format = "%.4f"
    f.write("#1.3\n")
    f.write(str(adata.shape[0]) + "\t" + str(adata.shape[1]) + "\t" + str(adata.obs.shape[1]) + "\t" + str(
        adata.var.shape[1]) + "\n")
    write_top_half_gct(f, adata.obs, adata.var, '', '')
    write_bottom_half_gct(f, adata.obs, pd.DataFrame(adata.X), '', data_float_format, '')

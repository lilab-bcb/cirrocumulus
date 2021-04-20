from urllib.parse import urlparse

import fsspec


def get_scheme(path):
    pr = urlparse(path)
    if len(pr.scheme) <= 1:  # for file paths: /foo/bar/test.h5ad or C:/foo/bar/test.h5ad
        return 'file'
    return pr.scheme


class FileSystemAdapter:

    def __init__(self):
        # scheme can be gs, file, etc.
        self.scheme_to_fs = {}

    def get_fs(self, path):
        scheme = get_scheme(path)
        fs = self.scheme_to_fs.get(scheme, None)
        if fs is not None:
            return fs
        fs = fsspec.filesystem(scheme)
        self.scheme_to_fs[scheme] = fs
        return fs

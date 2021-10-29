import os

from flask import current_app

from .envir import CIRRO_AUTH, CIRRO_DATABASE, CIRRO_MOUNT


def get_database():
    return current_app.config[CIRRO_DATABASE]


def get_auth():
    return current_app.config[CIRRO_AUTH]


remapped_urls = dict()
if os.environ.get(CIRRO_MOUNT) is not None:
    tokens = os.environ.get(CIRRO_MOUNT).split(',')
    for token in tokens:
        index = token.rfind(':')
        bucket = token[:index]
        local_path = token[index + 1:]
        remapped_urls[bucket] = local_path


def map_url(url):
    for remote_path in remapped_urls:
        before_replace = url
        url = url.replace(remote_path, remapped_urls[remote_path])
        if before_replace != url:  # only replace URL once
            return url
    return url

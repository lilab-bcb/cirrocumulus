import json
import logging
import os
from abc import abstractmethod, ABC

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse

from cirro.simple_data import SimpleData

logger = logging.getLogger("cirro")


class AbstractDataset(ABC):

    def __init__(self):
        super().__init__()
        self.cached_dataset_id = None
        self.cached_data = {}

    @abstractmethod
    def get_suffixes(self):
        pass

    @abstractmethod
    def to_pandas(self, file_system, path, columns):
        pass

    def schema(self, file_system, path):
        with file_system.open(path) as s:
            return json.load(s)

    @staticmethod
    def get_keys(keys, basis=None):
        if basis is not None:
            logger.info(keys)
            logger.info(basis)
            keys = keys + basis['coordinate_columns']
        return keys

    def read_summarized(self, file_system, path, obs_keys=[], var_keys=[], index=False, rename=False, dataset=None):
        result_df = pd.DataFrame()
        if index:
            df = self.to_pandas(file_system, path + '/index')
            for column in df:
                result_df[column] = df[column]

        for key in var_keys + obs_keys:
            df = self.to_pandas(file_system, path + '/' + key)
            if rename:
                for column in df:
                    result_df['{}_{}'.format(key, column)] = df[column]
            else:
                for column in df:
                    result_df[column] = df[column]
        return result_df

    def read_data_sparse(self, file_system, path, obs_keys=[], var_keys=[], basis=[], dataset=None):
        # path is path to index.json
        # if path ends with /data then X is stored as index, value pairs
        # if basis, read index to get bins, x, and y
        path = os.path.dirname(path)
        data_path = os.path.join(path, 'data')
        data = []
        row = []
        col = []
        X = None
        nobs = dataset['nObs']

        if len(var_keys) > 0:
            for i in range(len(var_keys)):
                key = var_keys[i]
                df = self.to_pandas(file_system, data_path + '/' + key)
                data.append(df['value'])
                row.append(df['index'])
                col.append(np.repeat(i, len(df)))

            data = np.concatenate(data)
            row = np.concatenate(row)
            col = np.concatenate(col)
            X = scipy.sparse.csr_matrix((data, (row, col)), shape=(nobs, len(var_keys)))
        obs = None
        dataset_id = dataset.id
        if self.cached_dataset_id != dataset_id:
            self.cached_dataset_id = dataset_id
            self.cached_data = {}

        for key in obs_keys:
            cache_key = str(dataset_id) + '-' + key
            cached_value = self.cached_data.get(cache_key)
            if cached_value is None:
                df = self.to_pandas(file_system, data_path + '/' + key, ['value'])  # ignore index in obs for now
                cached_value = df['value']
            if obs is None:
                obs = pd.DataFrame()
            obs[key] = cached_value
            self.cached_data[cache_key] = cached_value

        if basis is not None and len(basis) > 0:
            for b in basis:
                cache_key = str(dataset_id) + '-' + b['full_name']
                is_precomputed = b['precomputed']
                if is_precomputed:
                    columns_to_fetch = [b['full_name']]
                else:  # need coordinates
                    columns_to_fetch = b['coordinate_columns']
                cached_value = self.cached_data.get(cache_key)
                if cached_value is None:
                    basis_path = (data_path + '/') + (b['full_name' if is_precomputed else 'name'])
                    df = self.to_pandas(file_system, basis_path, columns_to_fetch)
                    cached_value = df
                    self.cached_data[cache_key] = cached_value
                if obs is None:
                    obs = pd.DataFrame()
                for c in columns_to_fetch:
                    obs[c] = cached_value[c]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

    def read(self, file_system, path, obs_keys=[], var_keys=[], basis=None, dataset=None):
        # path is path to index.json
        if path.endswith('json'):
            return self.read_data_sparse(file_system, path, obs_keys, var_keys, basis, dataset)
        keys = []
        if basis is not None and len(basis) > 0:
            for b in basis:
                keys += AbstractDataset.get_keys(obs_keys + var_keys, basis=b)
        if len(keys) == 0:
            return SimpleData(None, pd.DataFrame(), pd.Index([]))
        df = self.to_pandas(file_system, path, keys)
        X = df[var_keys].values
        if basis is not None and len(basis) > 0:
            for b in basis:
                obs_keys = obs_keys + b['coordinate_columns']
        obs = df[obs_keys]
        return SimpleData(X, obs, pd.DataFrame(index=pd.Index(var_keys)))

import gzip
import json
import unittest

import cirro.api


class TestEmbedding(unittest.TestCase):
    def test_one_feature(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=['col4'], view_column_names=['col5'],
            nbins=2)
        self.assertEqual(result['embedding']['col4'][0], 4.0)
        self.assertEqual(len(result['embedding']['col5']), 2)

    def test_two_feature(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=['col4', 'col5'], view_column_names=['col6'],
            nbins=2)
        self.assertEqual(result['embedding']['col4'][0], 4.0)

    def test_nos(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=['col4'],
            view_column_names=['col5'], nbins=None)
        self.assertEqual(len(result['embedding']['col5']), 5)

    def test_nos_no_features(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f,
            view_column_names=['col5'], nbins=None)
        self.assertEqual(len(result['embedding']['col5']), 5)

    def test_grouping(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, group_by_column_names=['col3'], feature_names=['col4'],
            view_column_names=['col5'], nbins=None)
        self.assertEqual(result['embedding']['col5'], [1, 2, 3, 4, 5])

    def test_grouping_one_feature(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, group_by_column_names=['col3'], feature_names=['col4'],
            view_column_names=['col5'], nbins=2)
        self.assertEqual(result['embedding']['col3_aa'][0], 2)

    def test_no_features(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=[], view_column_names=['col6'],
            nbins=2)
        self.assertEqual(result['embedding']['density'][0], 4)

    def test_gzip(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=[], view_column_names=['col6'],
            nbins=2)
        result = cirro.api.to_json(result)
        result = json.loads(gzip.decompress(result))
        self.assertEqual(result['embedding']['density'][0], 4)

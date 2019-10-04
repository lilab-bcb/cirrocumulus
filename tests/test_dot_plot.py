import unittest

import cirro.api


class TestDotPlot(unittest.TestCase):

    def test_one_feature(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=['col6'])
        self.assertEqual(result['stats']['var'][0]['mean'], 30.0)
        self.assertEqual(result['stats']['var'][0]['non_zero'], 1)

    def test_two_features(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, feature_names=['col6', 'col5'])
        self.assertEqual(result['stats']['var'][0]['mean'], 30.0)
        self.assertEqual(result['stats']['var'][0]['non_zero'], 1)
        self.assertEqual(result['stats']['var'][1]['mean'], 3.0)
        self.assertEqual(result['stats']['var'][1]['non_zero'], 1)

    def test_grouping_one_feature(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, group_by_column_names=['col3'], feature_names=['col6'])
        self.assertEqual(result['stats']['var'][0]['mean'][result['stats']['obs'].index('aa')], 100.0 / 3.0)

    def test_grouping_two_features(self):
        f = open('test.parquet', mode='rb')
        result = cirro.api.slice_from_file(f=f, group_by_column_names=['col3'], feature_names=['col6', 'col5'])
        self.assertEqual(result['stats']['var'][0]['mean'][result['stats']['obs'].index('aa')], 100.0 / 3.0)
        self.assertEqual(result['stats']['var'][1]['mean'][result['stats']['obs'].index('aa')], 10.0 / 3.0)

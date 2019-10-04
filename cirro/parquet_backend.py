import pyarrow
import pyarrow.parquet as pq


class ParquetBackend:

    def schema(self, file_system, path):
        with file_system.open(path) as f:
            pq_schema = pq.read_schema(f)
        names = pq_schema.names
        types = pq_schema.types
        result = {'version': '1'}

        var = []
        obs_cat = []
        str_type = pyarrow.lib.string()
        layout_to_dimensions = {}
        for i in range(len(names)):
            key = names[i]
            if key.startswith('X_') and (
                    key.endswith('_1') or key.endswith('_2') or key.endswith('_3')):  # assume layout (e.g. X_PCA_1)
                basis = key[0:len(key) - 2]
                dimension = int(key[len(key) - 1:])
                max_dim = layout_to_dimensions.get(basis, 0)
                if dimension > max_dim:
                    layout_to_dimensions[basis] = dimension
            elif types[i] == str_type:
                obs_cat.append(key)
            else:
                var.append(key)
        result['var'] = var
        result['obs'] = []
        result['obs_cat'] = obs_cat
        layouts = []
        for key in layout_to_dimensions:
            layouts.append({'name': key, 'dimensions': layout_to_dimensions[key]})
        result['layouts'] = layouts
        return result

    def get_df(self, file_system, path, keys, layout_key=None):
        if layout_key is not None:
            layout_name = layout_key['name']
            for i in range(layout_key['dimensions']):
                keys.append(layout_name + '_' + str(i + 1))
        with file_system.open(path) as f:
            table = pq.read_table(f, columns=keys)
        return table.to_pandas()

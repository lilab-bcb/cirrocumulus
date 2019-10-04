import pandas as pd

import pyarrow as pa
import pyarrow.parquet as pq

df = pd.DataFrame(data={'col1': ['a', 'b', 'c', 'd', 'e'], 'col2': ['a1', 'b1', 'c1', 'd1', 'e1'],
                        'col3': ['qq', 'aa', 'aa', 'bb', 'aa'], 'col4': [1, 4, 5, 6, 7], 'col5': [1, 2, 3, 4, 5],
                        'col6': [10, 20, 30, 40, 50]})
table = pa.Table.from_pandas(df)
pq.write_table(table, 'test.parquet')

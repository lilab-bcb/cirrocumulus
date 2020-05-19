import anndata
import pandas as pd


def read_star_fusion_file(input_csv: str):
    """Load data from a STAR Fusion TSV file

    Parameters
    -----------

    input_csv : `str`
        The TSV file, gzipped or not, containing the STAR Fusion results.

    Returns
    -------
    """
    df = pd.read_csv(input_csv, sep='\t')
    df = df.sort_values('JunctionReadCount', ascending=False)
    df['value'] = True
    # df['value'] = df['value'].astype('category')
    df = df.drop_duplicates(['Cell', '#FusionName'], keep='first')
    # m = df.pivot_table('value', 'Cell', '#FusionName', fill_value=0)
    m = df.pivot('Cell', '#FusionName', 'value')
    m.fillna(False, inplace=True)
    return anndata.AnnData(obs=m)

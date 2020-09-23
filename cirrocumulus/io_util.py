import os

import anndata
import pandas as pd


def get_cell_type_genes(cell_type):
    all_genes = []

    for cell_type_markers in cell_type['markers']:
        all_genes += cell_type_markers['genes']

    for i in range(len(all_genes)):
        gene = all_genes[i]
        last_char = gene[len(gene) - 1]
        if last_char == '+' or last_char == '-':
            all_genes[i] = gene[:len(gene) - 1]
    return all_genes


# remove genes in dict that are not in dataset
def filter_markers(adata, marker_dict):
    for category_name in marker_dict:
        markers = marker_dict[category_name]
        for key in markers:
            prior_size = len(markers[key])
            filtered_markers = set(markers[key]).intersection(adata.var.index)
            if len(filtered_markers) != prior_size:
                markers[key] = list(filtered_markers)
    return marker_dict


def get_markers(marker_paths):
    marker_dict = {}
    import json

    for marker_path in marker_paths:
        markers = {}
        with open(marker_path, 'rt') as f:
            marker_json = json.load(f)
            if 'title' in marker_json:
                marker_dict[marker_json['title']] = markers
                for cell_type in marker_json['cell_types']:
                    markers[cell_type['name']] = get_cell_type_genes(cell_type)
                    if 'subtypes' in cell_type:
                        for cell_sub_type in cell_type['subtypes']['cell_types']:
                            markers[cell_sub_type['name']] = get_cell_type_genes(cell_sub_type)
            else:
                key = os.path.splitext(os.path.basename(marker_path))[0]
                markers = marker_dict.get(key, {})
                marker_dict[key] = markers
                markers.update(marker_json)
    return marker_dict


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

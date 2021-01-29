import os

import pandas as pd

SPATIAL_HELP = 'Directory containing 10x visium spatial data (tissue_hires_image.png, scalefactors_json.json, ' \
               + 'and tissue_positions_list.csv) ' + 'or a directory containing `image.png`, `positions.image.csv` ' \
               + 'with headers barcode, x, and y, and optionally `diameter.image.txt` containing spot diameter'


def __add_visium(adata, spatial_directory):
    scale_factors_path = os.path.join(spatial_directory, 'scalefactors_json.json')
    tissue_hires_image_path = os.path.join(spatial_directory, 'tissue_hires_image.png')
    tissue_positions_list_path = os.path.join(spatial_directory, 'tissue_positions_list.csv')
    is_visium = True
    for path in [scale_factors_path, tissue_hires_image_path, tissue_positions_list_path]:
        if not os.path.exists(path):
            is_visium = False
            break
    if is_visium:
        import json
        with open(os.path.join(spatial_directory, 'scalefactors_json.json'), 'rt') as f:
            scalefactors = json.load(f)
            # {"spot_diameter_fullres": 89.49502418224989, "tissue_hires_scalef": 0.17011142,
            # "fiducial_diameter_fullres": 144.56888521748058, "tissue_lowres_scalef": 0.051033426}
        # barcode, in_tissue, array_row, array_col, pxl_col_in_fullres, pxl_row_in_fullres
        positions = pd.read_csv(tissue_positions_list_path, header=None)
        positions.columns = [
                'barcode',
                'in_tissue',
                'array_row',
                'array_col',
                'pxl_col_in_fullres',
                'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']
        positions = positions.reindex(adata.obs.index)
        spatial_coords = positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
        adata.obsm['tissue_hires'] = spatial_coords * scalefactors['tissue_hires_scalef']
        adata.uns['images'] = [dict(type='image', name='tissue_hires', image=tissue_hires_image_path,
            spot_diameter=scalefactors['spot_diameter_fullres'] * scalefactors['tissue_hires_scalef'])]
        return True
    else:
        return False


cirro_id_counter = 0


def cirro_id():
    global cirro_id_counter
    cirro_id_counter += 1
    return 'cirro-{}'.format(cirro_id_counter)


def unique_id():
    import uuid
    return str(uuid.uuid4())


def __add_generic_spatial(adata, spatial_directory):
    # positions.image.csv with barcode, x, y
    # image.png
    # diameter.image.txt (optional)
    image_extensions = set(['.png', '.jpeg', '.jpg'])
    meta_image_extensions = set(['.svg'])
    found = False
    for f in os.listdir(spatial_directory):
        name, ext = os.path.splitext(f)
        ext = ext.lower()
        if ext in image_extensions:
            positions_path = os.path.join(spatial_directory, 'positions.' + name + '.csv')
            if os.path.exists(positions_path):
                diameter_path = os.path.join(spatial_directory, 'diameter.' + name + '.txt')
                spot_diameter = None
                if os.path.exists(diameter_path):
                    with open(diameter_path, 'rt') as diameter_in:
                        spot_diameter = float(diameter_in.readline().strip())

                positions = pd.read_csv(positions_path, index_col='barcode')
                if not pd.api.types.is_string_dtype(positions.index):
                    positions.index = positions.index.astype(str)
                positions = positions.reindex(adata.obs.index)
                adata.obsm[name] = positions[['x', 'y']].values
                images = adata.uns.get('images', [])
                images.append(dict(type='image', name=name, image=os.path.join(spatial_directory, f),
                    spot_diameter=spot_diameter))
                adata.uns['images'] = images
                found = True
        elif ext in meta_image_extensions:
            svg_path = os.path.join(spatial_directory, f)
            if os.path.exists(svg_path):
                svg_path = os.path.abspath(svg_path)
            import xml.etree.ElementTree as ET
            import json

            tree = ET.parse(svg_path)
            attrs = tree.getroot().attrib
            if 'data-group' in attrs and 'data-selection' in attrs:
                found = True
                images = adata.uns.get('meta_images', [])
                selection = attrs['data-selection'].replace("'", "\"")
                selection = json.loads(selection)

                images.append(
                    dict(type='meta_image', name=name, image=svg_path,
                        attrs=dict(group=attrs['data-group'], selection=selection)))
                adata.uns['meta_images'] = images
    return found


def add_spatial(adata, spatial_directory):
    if not __add_visium(adata, spatial_directory):
        return __add_generic_spatial(adata, spatial_directory)
    else:
        return True


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


# remove genes that are not in dataset
def filter_markers(adata, markers):
    for result in markers:
        features = result['features']
        prior_size = len(features)
        filtered_markers = set(features).intersection(adata.var.index)
        if len(filtered_markers) != prior_size:
            result['features'] = list(filtered_markers)
    return markers


def get_markers(marker_paths):
    import json
    marker_results = []
    for marker_path in marker_paths:

        with open(marker_path, 'rt') as f:
            marker_json = json.load(f)
            if 'title' in marker_json:
                title = marker_json['title']
                for cell_type in marker_json['cell_types']:
                    name = cell_type['name']
                    features = get_cell_type_genes(cell_type)
                    marker_results.append(dict(category=title, name=name, features=features))
                    if 'subtypes' in cell_type:
                        for cell_sub_type in cell_type['subtypes']['cell_types']:
                            subtype_name = cell_sub_type['name']
                            subtype_features = get_cell_type_genes(cell_sub_type)
                            marker_results.append(dict(category=title, name=subtype_name, features=subtype_features))
            else:
                key = os.path.splitext(os.path.basename(marker_path))[0]
                for name in marker_json:
                    marker_results.append(dict(category=key, name=name, features=marker_json[name]))
    return marker_results


def read_star_fusion_file(input_csv: str):
    """Load data from a STAR Fusion TSV file

    Parameters
    -----------

    input_csv : `str`
        The TSV file, gzipped or not, containing the STAR Fusion results.

    Returns
    -------
    """
    import anndata
    df = pd.read_csv(input_csv, sep='\t')
    df = df.sort_values('JunctionReadCount', ascending=False)
    df['value'] = True
    # df['value'] = df['value'].astype('category')
    df = df.drop_duplicates(['Cell', '#FusionName'], keep='first')
    # m = df.pivot_table('value', 'Cell', '#FusionName', fill_value=0)
    m = df.pivot('Cell', '#FusionName', 'value')
    m.fillna(False, inplace=True)
    return anndata.AnnData(obs=m)

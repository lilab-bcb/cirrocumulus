import os.path
import argparse

import numpy as np
import anndata

import cirrocumulus.io_util
from cirrocumulus.anndata_dataset import read_adata
from cirrocumulus.util import get_fs


def concat_spatial(paths: list[str], output_dir: str, ncols: int = 2):
    datasets = []
    images = []
    spot_diameters = []
    common_obsm_keys = None
    for path in paths:
        filesystem = get_fs(path)
        adata = read_adata(path, filesystem, False)
        adata.var.index = adata.var.index.str.upper()
        dataset_name = os.path.basename(path)
        dataset_name = os.path.splitext(dataset_name)[0]
        if dataset_name.endswith("_filtered_feature_bc_matrix"):
            dataset_name = dataset_name[: -len("_filtered_feature_bc_matrix")]

        spatial_dir = os.path.join(os.path.dirname(path), "spatial")
        if os.path.exists(spatial_dir) and cirrocumulus.io_util.add_spatial(adata, spatial_dir):
            spatial = adata.uns["images"]
            spot_diameters.append(spatial[0]["spot_diameter"])
            import imageio.v3 as iio

            images.append(iio.imread(spatial[0]["image"]))
        adata.obs["dataset"] = dataset_name
        if common_obsm_keys is None:
            common_obsm_keys = set()
            common_obsm_keys.update(adata.obsm.keys())
        else:
            common_obsm_keys.intersection(adata.obsm.keys())
        adata.obs.index = adata.obs.index + "-" + dataset_name
        datasets.append(adata)

    ncols = min(ncols, len(datasets))
    nrows = len(datasets) // ncols
    max_width = 0
    max_height = 0
    for i in range(len(images)):
        image_shape = images[i].shape
        max_width = max(max_width, image_shape[1])
        max_height = max(max_height, image_shape[0])
    indices = []
    for i in range(nrows):
        for j in range(ncols):
            indices.append((i, j))

    do_concat_spatial = len(images) == len(datasets)
    common_obsm = dict()
    if len(common_obsm_keys) > 0:
        for key in common_obsm_keys:
            if key != "tissue_hires":
                coords_list = []
                for i in range(len(datasets)):
                    coords = datasets[i].obsm[key]
                    if coords.shape[1] in [2, 3]:
                        row_index, col_index = indices[i]

                        coords[:, 0] = np.interp(
                            coords[:, 0],
                            (coords[:, 0].min(), coords[:, 0].max()),
                            (col_index, col_index + 1),
                        )
                        coords[:, 1] = np.interp(
                            coords[:, 1],
                            (coords[:, 1].min(), coords[:, 1].max()),
                            (row_index, row_index + 1),
                        )
                        coords_list.append(coords)
                if len(coords_list) == len(datasets):
                    common_obsm[key] = np.concatenate(coords_list)
    if do_concat_spatial:
        combined_image = np.zeros(
            (max_height * nrows, max_width * ncols, images[0].shape[2]), dtype=images[0].dtype
        )
        combined_image[:] = images[0].max()
        adata_spatial_coords = []
        for i in range(len(datasets)):
            spatial_coords = datasets[i].obsm["tissue_hires"]  # x, y coords
            row_index, col_index = indices[i]
            if row_index > 0:
                spatial_coords[:, 1] = spatial_coords[:, 1] + max_height * row_index
            if col_index > 0:
                spatial_coords[:, 0] = spatial_coords[:, 0] + max_width * col_index
            adata_spatial_coords.append(spatial_coords)
            img = images[i]
            x_start = row_index * max_height
            y_start = col_index * max_width
            combined_image[x_start : x_start + img.shape[0], y_start : y_start + img.shape[1]] = img

    combined_adata = anndata.concat(datasets, join="outer")
    combined_adata.obs_names_make_unique()
    for key in common_obsm:
        combined_adata.obsm[key] = common_obsm[key]
    os.makedirs(output_dir, exist_ok=True)
    if do_concat_spatial:
        combined_adata.obsm["tissue_hires"] = np.concatenate(adata_spatial_coords)
        spatial_dir = os.path.join(output_dir, "spatial")
        os.makedirs(spatial_dir, exist_ok=True)
        combined_image_path = os.path.join(spatial_dir, "tissue_hires.png")
        iio.imwrite(combined_image_path, combined_image)
        _spot_diameter = spot_diameters[0]
        spot_diameters_differ = False
        for i in range(1, len(spot_diameters)):
            if spot_diameters[i] != _spot_diameter:
                spot_diameters_differ = True
                break
        if spot_diameters_differ:  # each spot has own diameter
            _spot_diameters = []
            for i in range(len(spot_diameters)):
                per_spot_diameter = np.zeros(datasets[i].shape[0])
                per_spot_diameter[:] = spot_diameters[i]
                _spot_diameters.append(per_spot_diameter)
            spot_diameter = np.concatenate(_spot_diameters)
            assert len(spot_diameter) == combined_adata.shape[0]
        else:
            spot_diameter = spot_diameters[0]
        combined_adata.uns["images"] = dict(
            type="image",
            name="tissue_hires",
            image=combined_image_path,
            spot_diameter=spot_diameter,
        )

    combined_adata.write(os.path.join(output_dir, "concat_data.h5ad"))


def create_parser(description=False):
    parser = argparse.ArgumentParser(
        description="Concatenate datasets in a grid layout. If all the datasets are spatial datasets, then tissue images are concatenated."
        if description
        else None
    )
    parser.add_argument(
        "dataset",
        help="Paths to dataset in h5ad, 10x h5, loom, or Seurat (rds) format. For spatial dataset, a 'spatial' directory must exist in the same directory as each dataset",
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory to write concatenated h5ad file and spatial directory if input datasets are spatial",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--cols",
        dest="cols",
        help="Number of columns in tiled layout",
        default=2,
        type=int,
        required=False,
    )
    return parser


def main():
    parser = create_parser(True)
    args = parser.parse_args()
    concat_spatial(args.dataset, output_dir=args.output, ncols=args.cols)


if __name__ == "__main__":
    main()

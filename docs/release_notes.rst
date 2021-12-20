Release Notes
-------------

- 1.1.25 `December 20, 2021`
    * Added option to compute differential expression results for all groups vs. rest
    * Updated tooltip display

- 1.1.24.post4 `December 2, 2021`
    * Removed --email flag from `cirro serve`
    * Fixed error storing large precomputed results

- 1.1.24.post3 `November 23, 2021`
    * Enable double-click to select on spatial images
    * Use spawn for running background jobs
    * Fixed bug in which dataset dialog remained open on save
    * Fixed error getting precomputed results

- 1.1.24.post2 `November 19, 2021`
    * Fixed error editing dataset

- 1.1.24.post1 `November 15, 2021`
    * Fixed error selecting all features in differential expression results
    * Added additional fields for new datasets

- 1.1.24 `November 10, 2021`
    * Changed default output format to zarr_ in `cirro prepare_data`
    * Fixed error retrieving job results stored in GridFS

- 1.1.23.post1 `November 9, 2021`
    * Fixed gzip encoding when writing to remote storage (e.g S3)

- 1.1.23 `November 3, 2021`
    * Style updates
    * Fixed error loading differential expression results

- 1.1.22.post10 `November 2, 2021`
    * Made left-side drawer collapsible
    * Fixed error saving cluster markers
    * Improved point size when saving embedding images

- 1.1.22.post9 `October 29, 2021`
    * Fixed error showing gallery labels

- 1.1.22.post8 `October 28, 2021`
    * Show zero cells instead of all cells when no cells pass filters

- 1.1.22.post7 `October 28, 2021`
    * Fixed default color scale for displaying features in sets
    * Changed single-click to double-click for selecting categorical values on primary embedding

- 1.1.22.post6 `October 25, 2021`
    * Preserve embeddings when multiple datasets passed to `cirro prepare_data`
    * Fixed error saving dataset views

- 1.1.22.post5 `October 22, 2021`
    * Separate color schemes for continuous observations, features, and modules
    * Explore gene modules stored in adata.uns['module']

- 1.1.22.post4 `October 19, 2021`
    * Added ability to select all members of a set
    * Fixed bug displaying renamed category labels on embedding
    * Added option to integrate Mixpanel_. Set environment variable CIRRO_MIXPANEL to your project token to track open dataset events

- 1.1.22.post3 `October 13, 2021`
    * Fixed bug that required two clicks to load differential expression results
    * Remove jobs from database when deleting datasets

- 1.1.22.post2 `October 12, 2021`
    * Fixed error when computing differential expression between two lasso'ed selections using `cirro serve`
    * Fixed display of "Sign In" button using `cirro serve`
    * Added ability to control the visibility of table columns in dataset chooser

- 1.1.22.post1 `October 8, 2021`
    * Save categorical legend scroll bar position when switching features

- 1.1.22 `October 7, 2021`
    * Added `--ontology` option to `cirro serve` and `cirro launch`
    * Fixed error saving colors from h5ad files in `cirro prepare`
    * Added option to search specific fields in dataset chooser
    * Save results in GridFS when using MongoDB without `--results` option

- 1.1.21 `October 5, 2021`
    * Removed `--backed` option in `cirro launch`
    * Added zarr_ format support in `cirro prepare_data`
    * Added `--results` option to `cirro serve` and `cirro launch`
    * Added ability to compute differential expression between all pairs of clusters
    * Select category by clicking cell on an embedding

- 1.1.20.post4 `September 17, 2021`
    * Style updates

- 1.1.20.post3 `September 10, 2021`
    * Fixed error computing differential expression in `cirro prepare_data` using `Pegasus`_/`Cumulus`_ when categories contain the : character

- 1.1.20.post2 `September 9, 2021`
    * `cirro prepare_data` can output directly to S3 or GCP bucket

- 1.1.20.post1 `September 8, 2021`
    * Save categorical colors to database
    * Save cluster positive and negative markers to database

- 1.1.20 `September 2, 2021`
    * Added static website hosting capabilities
    * Performance improvements

- 1.1.19.post1 `August 25, 2021`
    * Fixed error getting precomputed results

- 1.1.19 `August 25, 2021`
    * Show distributions for numerical cell metadata
    * Improve interactive differential expression performance

- 1.1.18 `August 16, 2021`
    * `cirro prepare_data` accepts multiple input datasets to better support multimodal data

- 1.1.17.post4 `August 13, 2021`
    * Fixed error in `cirro launch` that prevented h5ad files from loading

- 1.1.17.post3 `August 13, 2021`
    * Compute complete differential expression results in `cirro prepare_data` using `Scanpy`_ or `Pegasus`_/`Cumulus`_

- 1.1.17.post2 `August 11, 2021`
    * Preserve category order only when < 1000 categories
    * Synchronize 3-d gallery chart rotation with primary view
    * Added separate marker size for filtered points

- 1.1.17.post1 `June 29, 2021`
    * Fixed error computing differential expression results when using `cirro launch`

- 1.1.17 `June 28, 2021`
    * Added TileDB support

- 1.1.16.post3 `June 3, 2021`
    * Embedding chart performance improvements
    * Replace saved filters with links in order to save complete state

- 1.1.16.post2 `May 26, 2021`
    * Drag and drop chips to reorder
    * Handle thousands of categories in violin plot

- 1.1.16.post1 `April 28, 2021`
    * Plot selected cells on top of unselected cells in embedding chart

- 1.1.16 `April 28, 2021`
    * Enable selecting top markers by a field in ascending or descending order
    * Updated auto-display logic

- 1.1.15.post4 `April 22, 2021`
    * Changed `--header` flag in `cirro serve` to accept Markdown file

- 1.1.15.post3 `April 20, 2021`
    * Fixed opening files with drive names on Windows

- 1.1.15.post2 `April 13, 2021`
    * Fixed error adding new dataset on Google App Engine
    * Show cirrocumulus version

- 1.1.15.post1 `April 12, 2021`
    * Moved composition plot to separate tab
    * Added `--header` flag to `cirro serve` to customize application header
    * Auto-display cluster annotation by default

- 1.1.15 `April 6, 2021`
    * Added composition plot
    * Pass `--upload` flag to `cirro serve` to enable file uploads
    * Show plot tooltips in bottom bar
    * Export data from dot plot

- 1.1.14.post5 `March 30, 2021`
    * Fixed issue that distribution charts did not update when color scheme changed

- 1.1.14.post4 `March 26, 2021`
    * Fixed issue that primary embedding chart did not update when color scheme changed

- 1.1.14.post3 `March 9, 2021`
    * Added ability to customize footer in `cirro serve`

- 1.1.14.post2 `March 2, 2021`
    * Click and drag to resize primary embedding chart
    * Added landing page

- 1.1.14.post1 `February 24, 2021`
    * Fixed error performing interactive differential expression analysis using `cirro launch`
    * Sort gallery charts by first by feature and then by embedding

- 1.1.14 `February 23, 2021`
    * Added interactive differential expression analysis
    * To add to current selection, hold down the Ctrl or Command keys when using lasso or box select tools

- 1.1.13.post2 `February 10, 2021`
    * Added standardize option that scales each feature or categorical group from zero to one for distributions and results visualization
    * Added species to dataset import when using `cirro serve`
    * Added option to show/hide labels in embedding gallery
    * `cirro launch` now accepts `Seurat`_ objects

- 1.1.13.post1 `February 2, 2021`
    * Added sort functionality to full differential expression results visualization

- 1.1.13 `February 1, 2021`
    * Explore complete differential expression results generated by `Scanpy`_ or `Pegasus`_/`Cumulus`_
    * Added reverse option to color schemes

- 1.1.12 `January 20, 2021`
    * Added violin plots

- 1.1.11.post3 `December 14, 2020`
    * Include categorical labels and dot plot options in `Copy Link` URL

- 1.1.11.post2 `December 8, 2020`
    * Use `anndata.uns[field_colors]` if present for cell metadata default colors
    * Added ability to view features in saved sets
    * Use `reticulate` to convert Seurat objects to h5ad in `cirro prepare_data`

- 1.1.11.post1 `December 6, 2020`
    * Convert seurat_clusters cell metadata field in Seurat objects to categorical in `cirro prepare_data`

- 1.1.11 `December 4, 2020`
    * Automatically compute cluster markers when using `cirro prepare_data` without --group flag
    * Show categorical labels on gene/feature embedding plots
    * Updated code for reading Seurat objects in `cirro prepare_data`

- 1.1.10.post8 `November 24, 2020`
    * Fixed error in `cirro prepare_data` when saving cell metadata names containing spaces

- 1.1.10.post7 `November 23, 2020`
    * Plot higher values on top of lower values for continuous values in saved embedding image.
    * Improved performance computing markers using `cirro prepare_data` with --group flag

- 1.1.10.post6 `November 20, 2020`
    * Fixed bug that prevented genes in sets from being displayed in selection dot plot.

- 1.1.10.post5 `November 18, 2020`
    * Fixed error when computing markers using `cirro prepare_data` with --group flag
    * Added ability to enter dataset description in Markdown when using `cirro serve`

- 1.1.10.post4 `November 12, 2020`
    * Toggle between dot plot and heatmap

- 1.1.10.post3 `November 6, 2020`
    * Added option to change dot plot color scheme

- 1.1.10.post2 `October 30, 2020`
    * Fixed display of set names
    * Fixed bug updating selected dot plot when selection changes

- 1.1.10.post1 `October 28, 2020`
    * Create dot plots by grouping by more than one category
    * Search dataset names and descriptions when using `cirro serve`

- 1.1.10 `October 25, 2020`
    * Fixed error selecting more than one cell metadata field

- 1.1.9.post3 `October 21, 2020`
    * Fixed error on startup using `cirro launch`

- 1.1.9.post2 `October 20, 2020`
    * Fixed serving spatial images using `cirro serve`

- 1.1.9.post1 `October 13, 2020`
    * Fixed error reading old datasets generated with `cirro prepare_data`

- 1.1.9 `October 13, 2020`
    * Added user interface to create gene/feature sets

- 1.1.8.post5 `October 5, 2020`
    * Updated dataset chooser

- 1.1.8.post4 `October 2, 2020`
    * Added dataset descriptions

- 1.1.8.post3 `October 1, 2020`
    * Show labels in gallery
    * Updated dark mode

- 1.1.8.post2 `September 29, 2020`
    * Removed active list. Select a feature/category to view details and filter

- 1.1.8.post1 `September 25, 2020`
    * Shuffle plot order in embedding plot for categorical values
    * Fixed scrolling bug in active list

- 1.1.8 `September 24, 2020`
    * Added support for generic spatial data in addition to 10x visium
    * Made primary embedding chart responsive
    * Added option to set min and max of color scale
    * Updated gallery chart size
    * Updated `prepare_data` command
    * Changed dot plot default min to zero
    * Added option to dot plot mean and percent expressed scales

- 1.1.7.post3 `September 18, 2020`
    * Plot higher values on top of lower values in embedding plot

- 1.1.7.post2 `September 17, 2020`
    * Save state when toggling between datasets
    * Fixed bug in dot plot tooltips
    * Changed dot plot color scheme

- 1.1.7.post1 `September 2, 2020`
    * Fixed bug passing `markers` to `launch` command
    * `launch` command takes multiple datasets

- 1.1.7 `August 28, 2020`
    * Use median instead of mean for categorical label position on data
    * Fixed Safari embedding label shadow bug
    * Save pan and zoom values in link URL

- 1.1.6 `August 27, 2020`
    * Added option to set embedding label font size
    * Show shadow around embedding label

- 1.1.5.post3 `August 26, 2020`
    * Fixed embedding label and tooltip color in dark mode
    * Fixed embedding label font size

- 1.1.5.post2 `August 25, 2020`
    * Save additional chart options when copying link
    * Support multiple differential expression results produced by `Scanpy`_

- 1.1.5.post1 `August 24, 2020`
    * Fixed dot plot background color in dark mode

- 1.1.5 `August 24, 2020`
    * Allow dataset sharing within an email domain
    * Added additional 3-d chart options
    * Added dark theme
    * Added timeout to `serve` command
    * Support markers generated with `Pegasus`_

- 1.1.4 `August 17, 2020`
    * Added spatial support

- 1.1.3 `August 13, 2020`
    * Improved support for Google authentication in `serve` command

- 1.1.2.post2 `August 12, 2020`
    * Fixed bug in `prepare_data` for saving markers
    * Added gunicorn and pymongo to requirements

- 1.1.2.post1 `August 11, 2020`
    * Added pyarrow to requirements

- 1.1.2 `August 11, 2020`
    * Show separate dot plots for all cells and selected cells
    * Added support for renaming clusters
    * Added `prepare_data` command for generating cirrocumulus formatted files for viewing on the cloud
    * Added 'serve' command to serve multiple users and datasets

- 1.1.1 `July 24, 2020`
    * Load marker genes from h5ad or JSON file

- 1.1.0.post3 `July 17, 2020`
    * Fixed embedding hover formatting issue

- 1.1.0.post2 `July 16, 2020`
    * Fixed Safari bug that caused gallery images to be flipped
    * Improved performance loading local h5ad files

- 1.1.0.post1 `June 15, 2020`
    * Fixed bug that sometimes prevented dot plot from showing

- 1.1.0 `June 1, 2020`
    * Added support for STARFusion output
    * Include labels in saved image

- 1.0.1 `May 7, 2020`
    * Draw labels on embedding

- 1.0.0 `May 5, 2020`
    * Lasso and box selection

- 0.0.6.post2 `Mar 25, 2020`
    * Added tabs for navigation
    * Use pandas for serialization

- 0.0.6.post1 `Mar 20, 2020`
    * Improved chart performance

- 0.0.6 `Mar 19, 2020`
    * Gallery view

- 0.0.5 `Mar 19, 2020`
    * Export filters

- 0.0.4 `Jan 16, 2020`
    * Autorotate 3d embeddings

- 0.0.3.post2 `Jan 14, 2020`
    * Save local filters to file

- 0.0.3.post1 `Jan 9, 2020`
    * Support 3d embeddings

- 0.0.3 `Jan 9, 2020`
    * Added filters
    * Added launch command

- 0.0.2 `Nov 5, 2019`
    * Initial release


.. _Pegasus: http://pegasus.readthedocs.io/
.. _Scanpy: https://scanpy.readthedocs.io/
.. _Seurat: https://satijalab.org/seurat/
.. _Cumulus: https://cumulus.readthedocs.io/en/stable/cumulus.html
.. _zarr: https://zarr.readthedocs.io/
.. _Mixpanel: https://mixpanel.com/
Release Notes
-------------

- 1.1.9.post4 `October 25, 2020`
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
    * Changed dotplot default min to zero
    * Added option to dotplot mean and percent expressed scales

- 1.1.7.post3 `September 18, 2020`
    * Plot higher values on top of lower values in embedding plot

- 1.1.7.post2 `September 17, 2020`
    * Save state when toggling between datasets
    * Fixed bug in dotplot tooltips
    * Changed dotplot color scheme

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
    * Fixed dotplot background color in dark mode

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
    * Show separate dotplots for all cells and selected cells
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
    * Fixed bug that sometimes prevented dotplot from showing

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

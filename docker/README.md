Example deploying cirrocumulus server, MongoDB, and [mongo-express](https://github.com/mongo-express/mongo-express)
using docker-compose. Datasets are stored in ./cirro_data.

In a production environment, we recommend using [AWS DocumentDB](https://aws.amazon.com/documentdb/) or [MongoDB Atlas](https://www.mongodb.com/atlas/database) instead of hosting MongoDB (and mongo-express) yourself.


Example
-------


- Download and unzip example dataset:

    ```
    curl -o pbmc3k.zip https://raw.githubusercontent.com/lilab-bcb/cirrocumulus/master/docs/pbmc3k.zip
    unzip pbmc3k.zip
    ```

- Convert h5ad to zarr format:

    ```
    cirro prepare_data pbmc3k/pbmc3k.h5ad --markers pbmc3k/markers.json --out cirro_data/pbmc3k.zarr
    ```

- Create and start containers:

    ```
    docker-compose up -d
    ```

- Import the dataset (users can also import datasets via the cirro user interface):

    ```
    curl http://localhost:3000/api/dataset -X POST -F 'name=my_name' -F 'url=/cirro_data/pbmc3k.zarr' -F 'description=my_desc'  -F 'species=Mus musculus'
    ```

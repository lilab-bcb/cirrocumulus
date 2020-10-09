// general server stuff

import {API, getIdToken} from './actions';

export class RestServerApi {

    getUserPromise() {
        return fetch(API + '/user', {
            headers: {
                'Authorization': 'Bearer ' + getIdToken(),
                'Content-Type': 'application/json'
            }
        }).then(result => result.json());
    }


    getDatasetsPromise() {
        return fetch(API + '/datasets', {headers: {'Authorization': 'Bearer ' + getIdToken()}})
            .then(response => {
                return response.json();
            });
    }

    deleteDatasetPromise(datasetId) {
        return fetch(API + '/dataset',
            {
                body: JSON.stringify(
                    {id: datasetId}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            });
    }

    upsertDatasetPromise(datasetId, dataset) {
        let isEdit = datasetId != null;
        if (datasetId != null) {
            dataset.id = datasetId;
        }
        return fetch(API + '/dataset',
            {
                body: JSON.stringify(dataset),
                method: isEdit ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(importDatasetResponse => importDatasetResponse.json());
    }


    getDatasetPromise(datasetId) {
        return fetch(API + '/dataset?id=' + datasetId,
            {
                method: 'GET',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }


// categories and filters
    setCategoryNamePromise(data) {
        return fetch(API + '/category_name',
            {
                body: JSON.stringify(data),
                method: 'PUT',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            });
    }

    upsertDatasetFilterPromise(data, isUpdate) {
        return fetch(API + '/filter',
            {
                body: JSON.stringify(data),
                method: isUpdate ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(response => response.json());
    }

    upsertFeatureSet(data, isUpdate) {
        return fetch(API + '/feature_set',
            {
                body: JSON.stringify(data),
                method: isUpdate ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(response => response.json());
    }

    deleteFeatureSet(setId, datasetId) {
        return fetch(API + '/feature_set',
            {
                body: JSON.stringify({id: setId, ds_id: datasetId}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            });
    }


    exportDatasetFiltersPromise(datasetId) {
        return fetch(API + '/export_filters?id=' + datasetId, {
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        }).then(result => {
            if (!result.ok) {
                return null;
            }
            return result.text();
        });
    }

    deleteDatasetFilterPromise(filterId, datasetId) {
        return fetch(API + '/filter',
            {
                body: JSON.stringify({id: filterId, ds_id: datasetId}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            });
    }

    getDatasetFilterPromise(filterId, datasetId) {
        return fetch(API + '/filter?id=' + filterId + '&ds_id=' + datasetId,
            {
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(response => response.json());
    }

    getCategoryNamesPromise(datasetId) {
        return fetch(API + '/category_name?id=' + datasetId,
            {
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }

    getFiltersPromise(datasetId) {
        return fetch(API + '/filters?id=' + datasetId,
            {
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }

}

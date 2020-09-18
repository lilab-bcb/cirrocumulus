import {API, getIdToken} from './actions';


export function getSelectedIdsPromise(data) {
    return fetch(API + '/selected_ids',
        {
            body: JSON.stringify(data),
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        });
}

export function getFileUrl(datasetId, file) {
    return API + '/file?id=' + datasetId + '&file=' + file + '&access_token=' + getIdToken();
}

export function getSelectionPromise(data) {
    return fetch(API + '/selection',
        {
            body: JSON.stringify(data),
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        }).then(result => result.json());
}

export function getSchemaPromise(id) {
    return fetch(API + '/schema?id=' + id, {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
        return response.json();
    });
}

export function getEmbeddingPromise(data) {
    return fetch(API + '/embedding',
        {
            body: JSON.stringify(data),
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        }).then(r => r.json());
}

export function getGroupedStatsPromise(data) {
    return fetch(API + '/grouped_stats',
        {
            body: JSON.stringify(data),
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        }).then(r => r.json());

}


export function getStatsPromise(data) {
    return fetch(API + '/stats',
        {
            body: JSON.stringify(data),
            method: 'POST',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        }).then(r => r.json());
}

// categories and filters
export function setCategoryNamePromise(data) {
    return fetch(API + '/category_name',
        {
            body: JSON.stringify(data),
            method: 'PUT',
            headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
        });
}

export function upsertDatasetFilterPromise(data, isUpdate) {
    return fetch(API + '/filter',
        {
            body: JSON.stringify(data),
            method: isUpdate ? 'PUT' : 'POST',
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        }).then(response => response.json());
}

export function exportDatasetFiltersPromise(id) {
    return fetch(API + '/export_filters?id=' + id, {
        headers: {'Authorization': 'Bearer ' + getIdToken()},
    }).then(result => {
        if (!result.ok) {
            return null;
        }
        return result.text();
    });
}

export function deleteDatasetFilterPromise(filterId, datasetId) {
    return fetch(API + '/filter',
        {
            body: JSON.stringify({id: filterId, ds_id: datasetId}),
            method: 'DELETE',
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        });
}

export function getDatasetFilterPromise(filterId, datasetId) {
    return fetch(API + '/filter?id=' + filterId + '&ds_id=' + datasetId,
        {
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        }).then(response => response.json());
}

export function getCategoryNamesPromise(id) {
    return fetch(API + '/category_name?id=' + id,
        {
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        });
}

export function getFiltersPromise(id) {
    return fetch(API + '/filters?id=' + id,
        {
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        });
}





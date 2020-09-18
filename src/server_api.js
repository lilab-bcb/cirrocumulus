// general server stuff

import {API, getIdToken} from './actions';

export function getUserPromise() {
    return fetch(API + '/user', {
        headers: {
            'Authorization': 'Bearer ' + getIdToken(),
            'Content-Type': 'application/json'
        }
    }).then(result => result.json());
}


export function getServerPromise() {
    return fetch(API + '/server').then(result => result.json());
}


export function getDatasetsPromise() {
    return fetch(API + '/datasets', {headers: {'Authorization': 'Bearer ' + getIdToken()}})
        .then(response => {
            return response.json();
        });
}


export function deleteDatasetPromise(datasetId) {
    return fetch(API + '/dataset',
        {
            body: JSON.stringify(
                {id: datasetId}),
            method: 'DELETE',
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        });
}

export function upsertDatasetPromise(datasetId, url, name, readers) {
    let isEdit = datasetId != null;
    return fetch(API + '/dataset',
        {
            body: JSON.stringify(
                {
                    id: datasetId,
                    url: url,
                    name: name,
                    readers: readers
                }),
            method: isEdit ? 'PUT' : 'POST',
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        }).then(importDatasetResponse => importDatasetResponse.json());
}


export function getDatasetPromise(datasetId) {
    return fetch(API + '/dataset?id=' + datasetId,
        {
            method: 'GET',
            headers: {'Authorization': 'Bearer ' + getIdToken()},
        }).then(result => result.json());
}


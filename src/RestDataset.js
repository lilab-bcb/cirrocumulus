import {API, getIdToken} from './actions';

export class RestDataset {

    init(id) {
        this.id = id;
        return Promise.resolve();
    }

    getSelectedIdsPromise(data) {
        data.id = this.id;
        return fetch(API + '/selected_ids',
            {
                body: JSON.stringify(data),
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).result(result => result.json());
    }

    getFileUrl(file) {
        return API + '/file?id=' + this.id + '&file=' + file + '&access_token=' + getIdToken();
    }

    getDataPromise(data) {
        data.id = this.id;
        return fetch(API + '/data',
            {
                body: JSON.stringify(data),
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }

    getSchemaPromise() {
        return fetch(API + '/schema?id=' + this.id, {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        });
    }
}






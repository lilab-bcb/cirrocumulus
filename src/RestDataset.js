import {API, getIdToken} from './actions';
import {getPassingFilterIndices} from './dataset_filter';
import {cacheValues, computeDerivedStats} from './VectorUtil';

function reshapeDistributionResult(distribution) {
    const results = [];
    distribution.forEach(distributionResult => {
        const categories = distributionResult.categories;
        const dimension = distributionResult.name;
        for (let i = 0; i < distributionResult.values.length; i++) {
            results.push({
                dimension: dimension,
                name: categories[i],
                feature: distributionResult.values[i].name,
                mean: distributionResult.values[i].mean,
                percentExpressed: distributionResult.values[i].percentExpressed
            });
        }
    });
    return results;
}

export class RestDataset {

    /**
     *
     * @param id Dataset id
     * param url Dataset URL
     * @param local Whether stats should be computed locally or by the server
     * @returns {Promise<void>}
     */
    init(id, url, local = true) {
        this.id = id;
        this.local = local;
        return Promise.resolve();
    }

    getSelectedIdsPromise(q, cachedData) {
        q.id = this.id;
        if (this.local) {
            if (cachedData['index'] == null) {
                return this.getDataPromise({values: {dimensions: ['index']}}, cachedData).then(() => {
                    const selectedIndices = Array.from(getPassingFilterIndices(cachedData, q.filter));
                    const ids = [];
                    const index = cachedData['index'];
                    for (let i = 0, n = selectedIndices.length; i < n; i++) {
                        ids.push(index[selectedIndices[i]]);
                    }
                    return {ids: ids};
                });
            }
        } else {
            return fetch(API + '/selected_ids',
                {
                    body: JSON.stringify(q),
                    method: 'POST',
                    headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
                }).then(result => result.json());
        }
    }

    getFileUrl(file) {
        return API + '/file?id=' + this.id + '&file=' + file + '&access_token=' + getIdToken();
    }

    getJob(id, returnResults) {
        return fetch(API + '/job?id=' + id + '&ds_id=' + this.id + (!returnResults ? '&status=1' : ''), {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        });
    }

    getJobs(id) {
        return fetch(API + '/jobs?id=' + id, {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        });
    }

    deleteJob(id) {
        return fetch(API + '/job',
            {
                body: JSON.stringify(
                    {id: id}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            });
    }

    getDataPromise(data, cachedData) {
        data.id = this.id;
        let dataSend = data;
        const local = this.local;

        if (this.local) {
            dataSend = {};
            //  ['stats', 'groupedStats', 'embedding', 'selection', 'values'];
            if (data.embedding || data.values) {
                dataSend.id = this.id;
                dataSend.embedding = data.embedding;
                dataSend.values = data.values;
            }
        }
        let body = JSON.stringify(dataSend);
        let p = body !== '{}' ? fetch(API + '/data',
            {
                body: body,
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).then(r => r.json()).then(result => {
            // convert sparse to dense?

            // if (result.values) {
            //     for (let key in result.values) {
            //         let data = result.values[key];
            //         if (data.index) {  // sparse
            //             let values = new Float32Array(xxx);
            //             for (let i = 0, n = data.index.length; i < n; i++) {
            //                 values[data.index[i]] = data.value[i];
            //             }
            //             result.values[key] = values;
            //         }
            //     }
            // }
            cacheValues(result, cachedData);
            return result;
        }) : Promise.resolve({});
        return p.then(result => {
            if (local) {
                computeDerivedStats(result, data, cachedData);
            } else {
                if (result.distribution) {
                    result.distribution = reshapeDistributionResult(result.distribution);
                }
                if (result.selection && result.selection.distribution) {
                    result.selection.distribution = reshapeDistributionResult(result.selection.distribution);
                }
            }
            return result;
        });
    }

    getSchemaPromise() {
        return fetch(API + '/schema?id=' + this.id, {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        });
    }
}






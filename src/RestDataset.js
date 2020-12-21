import {API, getIdToken} from './actions';
import {cacheValues, computeDerivedStats} from './VectorUtil';

function reshapeDotPlotResult(dotplot) {
    const results = [];
    dotplot.forEach(dotPlotResult => {
        const categories = dotPlotResult.categories;
        const dimension = dotPlotResult.name;
        for (let i = 0; i < dotPlotResult.values.length; i++) {
            results.push({
                dimension: dimension,
                name: categories[i],
                feature: dotPlotResult.values[i].name,
                mean: dotPlotResult.values[i].mean,
                fractionExpressed: dotPlotResult.values[i].fractionExpressed
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

    getSelectedIdsPromise(data) {
        data.id = this.id;
        return fetch(API + '/selected_ids',
            {
                body: JSON.stringify(data),
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }

    getFileUrl(file) {
        return API + '/file?id=' + this.id + '&file=' + file + '&access_token=' + getIdToken();
    }

    getDataPromise(data, cachedData) {
        data.id = this.id;
        let dataSend = data;
        const local = this.local;
        let send = true;

        if (this.local) {
            dataSend = {};
            //  ['stats', 'groupedStats', 'embedding', 'selection', 'values'];
            if (data.embedding || data.values) {
                dataSend.id = this.id;
                dataSend.embedding = data.embedding;
                dataSend.values = data.values;
                send = true;
            }
        }
        let body = JSON.stringify(dataSend);
        let p = body !== '{}' ? fetch(API + '/data',
            {
                body: body,
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).then(r => r.json()).then(result => {
            // convert sparse to dense
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
                if (result.dotplot) {
                    result.dotplot = reshapeDotPlotResult(result.dotplot);
                }
                if (result.selection && result.selection.dotplot) {
                    result.selection.dotplot = reshapeDotPlotResult(result.selection.dotplot);
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






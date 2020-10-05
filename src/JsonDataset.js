import {isArray, isObject} from 'lodash';
import {getPassingFilterIndices} from './dataset_filter';
import {SlicedVector} from './SlicedVector';
import {Vector} from './Vector';
import {groupedStats, stats, valueCounts} from './VectorUtil';


function getVarNameType(key) {
    let index = key.indexOf('/');
    if (index === -1) {
        return {name: key, type: 'X'};
    } else {
        let key_type = key.substring(0, index);
        let name = key.substring(index + 1);
        return {name: name, type: key_type};
    }
}

export function getBasis(basis, nbins = null, agg = null, dimensions = 2, precomputed = false) {
    dimensions = parseInt(dimensions);
    let coordinate_columns = [];
    for (let i = 0; i < dimensions; i++) {
        coordinate_columns.push(basis + '_' + i + 1);
    }
    let full_name = basis + '_' + dimensions;
    if (nbins != null) {
        full_name = full_name + '_' + nbins + '_' + agg;
    }
    return {
        'name': basis, 'dimensions': dimensions, 'coordinate_columns': coordinate_columns, 'nbins': nbins,
        'agg': agg, 'full_name': full_name, 'precomputed': precomputed
    };
}

function splitDataFilter(data_filter) {
    let var_keys = new Set();
    let obs_keys = new Set();
    let basis_list = [];
    if (data_filter) {
        let user_filters = data_filter.filters || [];
        for (let i = 0; i < user_filters.length; i++) {
            let user_filter = user_filters[i];
            let key = user_filter[0];
            if (isObject(key)) {
                let basis = getBasis(key.basis, key.nbins, key.agg,
                    key.ndim || 2, key.precomputed);
                basis_list.push(basis);
            } else {

                const {name, type} = getVarNameType(key);

                user_filter[0] = name;
                if (type === 'X') {
                    var_keys.add(name);
                } else {
                    obs_keys.add(name);
                }
            }
        }
    }
    return {basis: basis_list, X: Array.from(var_keys), obs: Array.from(obs_keys)};
//    return list(var_keys), list(obs_keys), basis_list
}

function splitMeasures(measures) {
    let obsMeasures = [];
    let varMeasures = [];
    for (let i = 0; i < measures.length; i++) {
        const {name, type} = getVarNameType(measures[i]);
        if (type === 'X') {
            varMeasures.push(name);
        } else if (type === 'obs') {
            obsMeasures.push(name);
        } else {
            throw('Unknown key type ' + type);
        }
    }
    return {obsMeasures: obsMeasures, varMeasures: varMeasures};
}


function getStats(dimensions, obsMeasures, varMeasures) {
    const results = {};
    dimensions.forEach(v => {
        results[v.getName()] = valueCounts(v);
    });
    obsMeasures.forEach(v => {
        results[v.getName()] = stats(v);
    });
    varMeasures.forEach(v => {
        results[v.getName()] = stats(v);
    });
    return results;
}

export class JsonDataset {

    init(id, url) {
        this.id = id;
        this.key2vector = {}; // maps key to vector

        this.schema = null;
        if (!url.endsWith('.json')) {
            url = url + 'schema.json';
        }
        this.url = url;
        this.baseUrl = this.url.substring(0, this.url.lastIndexOf('/') + 1);
        // return new Promise((resolve, reject) => {
        //     fetch(url + '.idx.json').then(r => r.json()).then(result => {
        //         this.key2bytes = result.index;
        //     }).then(() => {
        //         fetch(url, this.getByteRange('schema')).then(response => {
        //             return response.json();
        //         }).then(result => {
        //             this.schema = result["schema"];
        //             resolve();
        //         });
        //     });
        // });
        return new Promise((resolve, reject) => {
            fetch(url).then(r => r.json()).then(result => {
                this.schema = result;
                resolve();
            });
        });
    }

    // getByteRange(key) {
    //     let range = this.key2bytes[key];
    //     if (!range) {
    //         throw key + ' not found';
    //     }
    //     return {headers: {'Range': 'bytes=' + range[0] + '-' + range[1]}};
    // }

    _fetch(key) {
        return fetch(this.baseUrl + key + '.json').then(r => r.json());
    }

    fetchData(keys) {
        let promises = [];
        keys.forEach(key => {
            if (this.key2vector[key] == null && key !== '__count') {
                let p = this._fetch(key).then(data => {
                    if (isArray(data)) {
                        this.key2vector[key] = new Vector(key, data);
                    } else {
                        // sparse
                        if (data.index) {
                            let values = new Float32Array(this.schema.shape[0]);
                            for (let i = 0, n = data.index.length; i < n; i++) {
                                values[data.index[i]] = data.value[i];
                            }
                            this.key2vector[key] = new Vector(key, values);
                        } else {
                            this.key2vector[key] = new Vector(key, data);
                        }
                    }

                });
                promises.push(p);
            }
        });
        return new Promise(resolve => {
            Promise.all(promises).then(() => resolve());
        });

    }

    getSelectedIdsPromise(q) {
        const dataFilter = q.filter;
        const {basis, X, obs} = splitDataFilter(dataFilter);
        let keys = [];
        basis.forEach(embedding => {
            keys.push(embedding.name);
        });
        keys.push('index');

        keys = keys.concat(X).concat(obs);
        return new Promise(resolve => {
            this.fetchData(keys).then(() => {
                let indices = getPassingFilterIndices(this.key2vector, dataFilter);
                let idVector = this.getVector('index', indices);
                let ids = [];
                for (let i = 0, n = idVector.size(); i < n; i++) {
                    ids.push(idVector.get(i));
                }
                resolve({ids: ids});
            });

        });
    }

    getDataPromise(q) {
        let dimensions = [];
        let measures = [];
        let queryKeys = ['stats', 'groupedStats', 'embedding', 'selection'];
        const results = {};
        queryKeys.forEach(key => {
            if (key in q) {
                let obj = q[key];
                if (!isArray(obj)) {
                    obj = [obj];
                }
                obj.forEach(value => {
                    if (value.dimensions) {
                        dimensions = dimensions.concat(value.dimensions);
                    }
                    if (value.measures) {
                        measures = measures.concat(value.measures);
                    }
                });

            }
        });

        const {obsMeasures, varMeasures} = splitMeasures(measures);
        let basisKeys = new Set();

        if (q.selection) { // get any embeddings
            const dataFilter = q.selection.filter;
            const {basis, X, obs} = splitDataFilter(dataFilter);

            dimensions = dimensions.concat(obs);
            measures = measures.concat(X);
            const embeddings = q.selection.embeddings || [];
            let mappedEmbeddings = [];
            embeddings.forEach(embedding => {
                let basis = getBasis(embedding.basis, embedding.nbins, embedding.agg,
                    embedding.ndim || 2, embedding.precomputed);
                basisKeys.add(basis.name);
                mappedEmbeddings.push(basis);
            });
            q.selection.embeddings = mappedEmbeddings;
            basis.forEach(embedding => {
                basisKeys.add(embedding.name);
            });

        }
        if (q.embedding) {
            q.embedding.forEach(embedding => {
                let basis = getBasis(embedding.basis, embedding.nbins, embedding.agg,
                    embedding.ndim || 2, embedding.precomputed);
                basisKeys.add(basis.name);
                embedding.basis = basis;
            });

        }


        return new Promise(resolve => {

            this.fetchData(dimensions.concat(obsMeasures).concat(varMeasures).concat(Array.from(basisKeys))).then(() => {

                if (q.embedding) {
                    results.embedding = [];
                    q.embedding.forEach(embedding => {
                        let dimensions = embedding.dimensions || [];
                        let measures = embedding.measures || [];
                        const {obsMeasures, varMeasures} = splitMeasures(measures);
                        let values = {};
                        let coordinates = this.getVector(embedding.basis.name).asArray();
                        dimensions.concat(obsMeasures).concat(varMeasures).forEach(key => {
                            if (key === '__count') {
                                values[key] = new Int8Array(coordinates[Object.keys(coordinates)[0]].length);
                                values[key].fill(1);
                            } else {
                                values[key] = this.getVector(key).asArray();
                            }
                        });
                        results.embedding.push({coordinates: coordinates, values: values});
                    });
                }
                if (q.stats) {
                    let dimensions = q.stats.dimensions || [];
                    let measures = q.stats.measures || [];
                    const {obsMeasures, varMeasures} = splitMeasures(measures);
                    results.summary = getStats(this.getVectors(dimensions), this.getVectors(obsMeasures), this.getVectors(varMeasures));
                }

                if (q.groupedStats) {
                    let dimensions = q.groupedStats.dimensions || [];
                    let measures = q.groupedStats.measures || [];
                    if (dimensions.length > 0 && measures.length > 0) {
                        results.dotplot = groupedStats(this.getVectors(dimensions), this.getVectors(measures));
                    }

                }

                if (q.selection) {
                    let dimensions = q.selection.dimensions || [];
                    let measures = q.selection.measures || [];
                    const embeddings = q.selection.embeddings || [];
                    const {obsMeasures, varMeasures} = splitMeasures(measures);
                    results.selection = {};

                    let indices = getPassingFilterIndices(this.key2vector, q.selection.filter);

                    if (embeddings.length > 0) {
                        results.selection.coordinates = {};
                        embeddings.forEach(embedding => {
                            results.selection.coordinates[embedding.full_name] = {'indices_or_bins': indices};
                        });
                    }

                    if (dimensions.length > 0 && varMeasures.length > 0) {
                        results.selection.dotplot = groupedStats(this.getVectors(dimensions, indices), this.getVectors(measures, indices));
                    }

                    results.selection.count = indices.length;
                    results.selection.summary = getStats(this.getVectors(dimensions, indices), this.getVectors(obsMeasures, indices), this.getVectors(varMeasures, indices));
                }
                resolve(results);
            });
        });
    }

    getSchemaPromise() {
        return Promise.resolve(this.schema);
    }

    getVector(key, indices = null) {
        let v = this.key2vector[key];
        return indices == null ? v : new SlicedVector(v, indices);
    }

    getFileUrl(file) {
        return this.baseUrl + file;
    }

    getVectors(keys, indices = null) {
        let result = [];
        keys.forEach(key => {
            let v = this.key2vector[key];
            if (indices != null) {
                v = new SlicedVector(v, indices);
            }
            result.push(v);
        });
        return result;
    }


}





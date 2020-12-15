import {isObject} from 'lodash';
import natsort from 'natsort';

import {getPassingFilterIndices} from './dataset_filter';
import {SlicedVector} from './SlicedVector';
import {Vector} from './Vector';

export function getVarNameType(key) {
    let index = key.indexOf('/');
    if (index === -1) {
        return {name: key, type: 'X'};
    } else {
        let key_type = key.substring(0, index);
        let name = key.substring(index + 1);
        return {name: name, type: key_type};
    }
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

function getVectors(cachedData, names, indices = null) {
    let vectors = [];
    names.forEach(name => {
        let array = cachedData[name];
        if (array == null) {
            throw name + " not found";
        }
        let v = new Vector(name, array);
        if (indices != null) {
            v = new SlicedVector(v, indices);
        }
        vectors.push(v);
    });
    return vectors;
}

export function cacheValues(result, cachedData) {
    if (result.embeddings != null) {
        result.embeddings.forEach(embedding => {
            if (embedding.coordinates && Object.keys(embedding.coordinates).length > 0) {
                cachedData[embedding.name] = embedding.coordinates;
            }
            if (embedding.values) { // binned values
                for (let feature in embedding.values) {
                    cachedData[feature + '_' + embedding.name] = embedding.values[feature];
                }
            }
        });
    }
    if (result.values) {
        for (let feature in result.values) {
            cachedData[feature] = result.values[feature];
        }
    }
}


export function splitMeasures(measures) {
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

export function splitDataFilter(data_filter) {
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


export function computeDerivedStats(result, q, cachedData) {

    if (q.stats) {
        const dimensions = q.stats.dimensions || [];
        const measures = q.stats.measures || [];
        const {obsMeasures, varMeasures} = splitMeasures(measures);

        result.summary = getStats(getVectors(cachedData, dimensions), getVectors(cachedData, obsMeasures), getVectors(cachedData, varMeasures));
    }

    if (q.groupedStats) {
        const dimensions = q.groupedStats.dimensions || []; // array of arrays
        const measures = q.groupedStats.measures || [];

        if (dimensions.length > 0 && measures.length > 0) {
            result.dotplot = groupedStats(getVectors(cachedData, dimensions[0]), getVectors(cachedData, measures));
        }
    }

    if (q.selection) {
        const dimensions = q.selection.dimensions || [];
        const measures = q.selection.measures || [];
        const embeddings = q.selection.embeddings || [];
        const {obsMeasures, varMeasures} = splitMeasures(measures);
        result.selection = {};

        let selectedIndices = getPassingFilterIndices(cachedData, q.selection.filter);
        if (embeddings.length > 0) {
            result.selection.coordinates = {};
            embeddings.forEach(embedding => {
                let basis = getBasis(embedding.basis, embedding.nbins,
                    embedding.agg, embedding.ndim || 2, embedding.precomputed);
                result.selection.coordinates[basis.full_name] = {'indices_or_bins': selectedIndices};
            });
        }

        if (dimensions.length > 0 && varMeasures.length > 0) {
            result.selection.dotplot = groupedStats(getVectors(cachedData, dimensions, selectedIndices), getVectors(cachedData, measures, selectedIndices));
        }
        result.selection.count = selectedIndices.length;
        result.selection.summary = getStats(getVectors(cachedData, dimensions, selectedIndices), getVectors(cachedData, obsMeasures, selectedIndices),
            getVectors(cachedData, varMeasures, selectedIndices));
    }
}


export function valueCounts(v) {
    let valueToCount = {};
    for (let i = 0, size = v.size(); i < size; i++) {
        let value = v.get(i);
        let count = valueToCount[value] || 0;
        valueToCount[value] = count + 1;
    }
    let keys = Object.keys(valueToCount);
    keys.sort(natsort());
    let counts = keys.map(key => valueToCount[key]);
    return {categories: keys, counts: counts};
}


export function groupedStats(dimensions, varMeasures) {
    const categoryToIndices = {};
    const ndim = dimensions.length;
    const size = dimensions[0].size();
    const tmp = [];
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < ndim; j++) {
            tmp[j] = dimensions[j].get(i);
        }
        const key = tmp.join(', ');
        let indices = categoryToIndices[key];
        if (indices === undefined) {
            indices = [];
            categoryToIndices[key] = indices;
        }
        indices.push(i);
    }
    let categories = Object.keys(categoryToIndices);
    categories.sort(natsort());
    let dimensionName = [];
    for (let j = 0; j < ndim; j++) {
        dimensionName[j] = dimensions[j].getName();
    }
    dimensionName = dimensionName.join('-');
    // each entry {dimension:dimensionName, name:category, feature:'', mean:0, fractionExpressed:xx}
    let result = [];
    categories.forEach(category => {
        const indices = categoryToIndices[category];
        varMeasures.forEach((v) => {

            // let otherIndices = [];
            // categories.forEach(otherCategory => {
            //     if (category !== otherCategory) {
            //         otherIndices = otherIndices.concat(categoryToIndices[otherCategory]);
            //     }
            // });
            // const restVector = new SlicedVector(v, otherIndices);
            // const restStats = stats(restVector);
            // results.values[index].restMean.push(restStats.mean);
            // results.values[index].restFractionExpressed.push(restStats.numExpressed / restVector.size());
            const categoryVector = new SlicedVector(v, indices);
            const categoryStats = stats(categoryVector);
            let y = null;
            let storeValues = false;
            if (storeValues) {
                y = new Float32Array(categoryVector.size());
                for (let i = 0, size = categoryVector.size(); i < size; i++) {
                    y[i] = categoryVector.get(i);
                }
            }
            // results.values[index].mean.push(categoryStats.mean);
            // results.values[index].fractionExpressed.push(categoryStats.numExpressed / categoryVector.size());
            const entry = {
                dimension: dimensionName,
                name: category,
                feature: v.getName(),
                mean: categoryStats.mean,
                fractionExpressed: categoryStats.numExpressed / categoryVector.size(),
                y: y
                //  plotly stuff
                // type: 'violin',
                // points: 'none',
                // spanmode: 'hard',
                // box: {
                //     line: {color: 'black'},
                //     visible: true
                // }
            };

            result.push(entry);
        });
    });
    return result;

}

function variance(v, mean) {
    const size = v.size();
    if (size <= 1) {
        return NaN;
    }
    let ss = 0;
    for (let j = 0; j < size; j++) {
        let x = v.get(j);
        let diff = x - mean;
        diff = diff * diff;
        ss += diff;
    }
    let n = size - 1;
    if (n < 1) {
        n = 1;
    }
    return ss / n;
}

export function stats(v) {
    let min = Number.MAX_VALUE;
    let max = -Number.MAX_VALUE;
    let mean = 0;
    let numExpressed = 0;
    if (v.size() === 0) {
        return {
            'min': Number.NaN, 'max': Number.NaN, 'sum': Number.NaN, 'mean': Number.NaN, 'numExpressed': Number.NaN,
        };
    }
    for (let i = 0, size = v.size(); i < size; i++) {
        const value = v.get(i);
        min = value < min ? value : min;
        max = value > max ? value : max;
        mean += value;
        if (value !== 0) {
            numExpressed++;
        }
    }
    let sum = mean;
    mean = mean / v.size();
    return {
        'min': min, 'max': max, 'sum': sum, 'mean': mean, 'numExpressed': numExpressed,
    };
}




import {isObject} from 'lodash';

import {getPassingFilterIndices} from './dataset_filter';
import {SlicedVector} from './SlicedVector';
import {FEATURE_TYPE, NATSORT} from './util';
import {Vector} from './Vector';

export function getVarNameType(key) {
    let index = key.indexOf('/');
    if (index === -1) {
        return {name: key, type: FEATURE_TYPE.X}; // default
    } else {
        const type = key.substring(0, index);
        const name = key.substring(index + 1);
        return {name: name, type: type};
    }
}

function getStats(dimensionVectors, measureVectors) {
    const results = {};
    dimensionVectors.forEach(v => {
        results[v.getName()] = valueCounts(v);
    });
    measureVectors.forEach(v => {
        results[v.getName()] = continuousVectorStats(v);
    });

    return results;
}


export function getVectors(cachedData, names, indices = null) {
    let vectors = [];
    names.forEach(name => {
        let array = cachedData[name];
        if (array == null) {
            throw new Error(name + " not found");
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


export function getTypeToMeasures(measures) {
    let typeToMeasures = {X: [], obs: []};
    for (let i = 0; i < measures.length; i++) {
        const {name, type} = getVarNameType(measures[i]);
        let typeMeasures = typeToMeasures[type];
        if (typeMeasures === undefined) {
            typeMeasures = [];
            typeToMeasures[type] = typeMeasures;
        }
        typeMeasures.push(name);

    }
    return typeToMeasures;
}

export function getBasis(basis, dimensions = 2, mode = null) {
    dimensions = parseInt(dimensions);
    let coordinate_columns = [];
    for (let i = 0; i < dimensions; i++) {
        coordinate_columns.push(basis + '_' + i + 1);
    }
    let full_name = basis + '_' + dimensions;
    if (mode != null) {
        full_name = full_name + '_' + mode;
    }
    return {
        name: basis,
        dimensions: dimensions,
        coordinate_columns: coordinate_columns,
        mode: mode,
        full_name: full_name
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
                basis_list.push(key);
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
        const typeToMeasures = getTypeToMeasures(measures);

        let measureVectors = [];
        Object.values(typeToMeasures).forEach(fields => {
            measureVectors = measureVectors.concat(getVectors(cachedData, fields));
        });

        result.summary = getStats(getVectors(cachedData, dimensions), measureVectors);
    }

    if (q.groupedStats) {
        const dimensions = q.groupedStats.dimensions || []; // array of arrays to support multiple groupings in the future
        const measures = q.groupedStats.measures || [];
        const typeToMeasures = getTypeToMeasures(measures);
        if (dimensions.length > 0) {
            // TODO, currently we only handle dimensions[0]
            const groupDimensionInfo = groupDimensions(getVectors(cachedData, dimensions[0]));
            const distribution = {};
            for (const key in typeToMeasures) {
                const measures = typeToMeasures[key];
                if (measures.length > 0) {
                    distribution[key] = groupedStats(groupDimensionInfo, getVectors(cachedData, measures));
                }
            }
            result.distribution = distribution;
        }
    }

    if (q.selection) {
        const dimensions = q.selection.dimensions || [];
        const measures = q.selection.measures || [];
        const typeToMeasures = getTypeToMeasures(measures);
        result.selection = {};
        result.selection.indices = getPassingFilterIndices(cachedData, q.selection.filter);
        const selectedIndices = Array.from(result.selection.indices);
        const dimensionVectors = getVectors(cachedData, dimensions, selectedIndices);
        if (dimensions.length > 0) {
            const groupDimensionInfo = groupDimensions(dimensionVectors);
            const distribution = {};
            for (const key in typeToMeasures) {
                const measures = typeToMeasures[key];
                if (measures.length > 0) {
                    distribution[key] = groupedStats(groupDimensionInfo, getVectors(cachedData, measures, selectedIndices));
                }
            }
            result.selection.distribution = distribution;
        }
        let measureVectors = [];
        Object.values(typeToMeasures).forEach(fields => {
            measureVectors = measureVectors.concat(getVectors(cachedData, fields, selectedIndices));
        });
        result.selection.summary = getStats(dimensionVectors, measureVectors);
    }
}


export function valueCounts(v) {
    let valueToCount = new Map();
    for (let i = 0, size = v.size(); i < size; i++) {
        const key = v.get(i);
        let count = valueToCount.get(key) || 0;
        valueToCount.set(key, count + 1);
    }
    const keys = Array.from(valueToCount.keys());
    keys.sort(NATSORT);
    const counts = keys.map(key => valueToCount.get(key));
    return {categories: keys, counts: counts};
}


export function groupDimensions(dimensions) {
    const keyToIndices = {};
    const ndim = dimensions.length;
    const size = dimensions[0].size();
    const tmp = [];
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < ndim; j++) {
            tmp[j] = dimensions[j].get(i);
        }
        const key = tmp.join(', ');
        let value = keyToIndices[key];
        if (value === undefined) {
            value = {key: tmp.slice(), indices: []};
            keyToIndices[key] = value;
        }
        value.indices.push(i);
    }
    let categories = Object.keys(keyToIndices);
    const dimensionNames = [];
    for (let j = 0; j < ndim; j++) {
        dimensionNames[j] = dimensions[j].getName();
    }
    const dimensionName = dimensionNames.join('-');
    return {categories, keyToIndices, dimensionName, dimensionNames};
}

export function groupedStats(groupDimensionInfo, vectors) {
    let result = [];
    groupDimensionInfo.categories.forEach(category => {
        const value = groupDimensionInfo.keyToIndices[category];
        vectors.forEach((v) => {
            const categoryVector = new SlicedVector(v, value.indices);
            const categoryStats = continuousVectorStats(categoryVector);
            const entry = {
                dimension: groupDimensionInfo.dimensionName,
                dimensions: groupDimensionInfo.dimensionNames,
                categories: value.key,
                feature: v.getName(),
                mean: categoryStats.mean,
                n: categoryVector.size(),
                percentExpressed: 100 * (categoryStats.numExpressed / categoryVector.size()),
                vector: categoryVector
            };
            result.push(entry);
        });
    });
    return result;

}

export function toArray(v) {
    const result = new Array(v.size());
    for (let i = 0, n = v.size(); i < n; i++) {
        result[i] = v.get(i);
    }
    return result;
}

export function toFloatArray(v) {
    const result = new Float32Array(v.size());
    for (let i = 0, n = v.size(); i < n; i++) {
        result[i] = v.get(i);
    }
    return result;
}

export function variance(v, mean) {
    const size = v.size();
    if (size <= 1) {
        return Number.NaN;
    }
    let ss = 0;
    let count = 0;
    for (let j = 0; j < size; j++) {
        let x = v.get(j);
        if (!Number.isNaN(x)) {
            let diff = x - mean;
            diff = diff * diff;
            ss += diff;
            count++;
        }
    }
    let n = count - 1;
    if (n < 1) {
        n = 1;
    }
    return ss / n;
}


export function continuousVectorStats(v) {

    if (v.size() === 0) {
        return {
            min: Number.NaN,
            max: Number.NaN,
            sum: Number.NaN,
            mean: Number.NaN,
            numExpressed: Number.NaN,
            logSum: Number.NaN
        };
    }
    let min = Number.MAX_VALUE;
    let max = -Number.MAX_VALUE;
    let sum = 0;
    let numExpressed = 0;
    let logSum = 0;
    let count = 0;
    for (let i = 0, size = v.size(); i < size; i++) {
        const value = v.get(i);
        if (!Number.isNaN(value)) {
            count++;
            min = value < min ? value : min;
            max = value > max ? value : max;
            sum += value;
            logSum += Math.expm1(value);
            if (value !== 0) {
                numExpressed++;
            }
        }
    }

    const mean = sum / count;
    return {
        min: min, max: max, sum: sum, mean: mean, numExpressed: numExpressed, n: v.size(), logSum: logSum
    };
}




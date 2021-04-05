import {isObject} from 'lodash';

import {getBasis, getVarNameType} from './VectorUtil';

function combine(a, b, op) {
    return op === 'or' ? new Set([...a, ...b]) : new Set(
        [...a].filter(x => b.has(x)));
}

function getIndices(array, f) {
    let result = new Set();
    for (let i = 0, size = array.length; i < size; i++) {
        if (f(array[i])) {
            result.add(i);
        }
    }
    return result;
}

export function createFilterFunction(filter) {
    const op = filter[1];
    let value = filter[2];
    let applyFunction;
    if (op === 'in') {
        value = new Set(value);
        applyFunction = (d) => value.has(d);
    } else if (op === '>') {
        applyFunction = (d) => d > value;
    } else if (op === '=') {
        applyFunction = (d) => d === value;
    } else if (op === '<') {
        applyFunction = (d) => d < value;
    } else if (op === '!=') {
        applyFunction = (d) => d !== value;
    } else if (op === '>=') {
        applyFunction = (d) => d >= value;
    } else if (op === '<=') {
        applyFunction = (d) => d <= value;
    } else {
        throw new Error('Unknown filter: ' + op);
    }
    return applyFunction;
}

export function getPassingFilterIndices(cachedData, data_filter) {
    let passingIndices = null;
    if (data_filter) {
        let user_filters = data_filter.filters || [];
        let combine_filters = data_filter.combine || 'and';
        for (let i = 0; i < user_filters.length; i++) {
            let filterObject = user_filters[i];
            let filterField = filterObject[0];
            let filterValue = filterObject[2];
            let keep = null;

            if (isObject(filterField)) { // selection box or lasso
                let selected_points_basis = getBasis(filterField['basis'], filterField.nbins,
                    filterField.agg, filterField.ndim || 2, filterField.precomputed);
                let coordinate_columns = selected_points_basis.coordinate_columns;
                if (filterValue.indices) { // Set of passing indices
                    let field = selected_points_basis['nbins'] ? selected_points_basis['full_name'] : 'index';
                    if (field == '__index') {
                        keep = filterValue.indices;
                    } else { // binning
                        throw new Error('Not implemented');
                    }
                    // keep = getIndices(cachedData[field], (val) => p.has(val));
                } else {
                    let selection_keep;
                    let path = filterValue.path;
                    for (let j = 0; j < path.length; j++) {
                        let p = path[j];
                        let xKeep = getIndices(cachedData[coordinate_columns[0]], (val) => val >= p.x && val <= p.x + p.width);
                        let yKeep = getIndices(cachedData[coordinate_columns[1]], (val) => val >= p.y && val <= p.y + p.height);
                        selection_keep = combine(xKeep, yKeep, 'and');
                        if (p.z) {  // 3d
                            let zKeep = getIndices(cachedData[coordinate_columns[2]], (val) => val >= p.z && val <= p.z + p.depth);
                            selection_keep = combine(selection_keep, zKeep, 'and');
                        }
                    }
                    keep = keep ? combine(selection_keep, keep, combine_filters) : selection_keep;
                }
            } else {
                if (filterField === '__index') {
                    keep = new Set(filterValue); // [__index, in, indices]
                } else {
                    const nameType = getVarNameType(filterField);
                    let series = cachedData[nameType.name];
                    keep = getIndices(series, createFilterFunction(filterObject));
                }
            }


            if (passingIndices) {
                passingIndices = combine(passingIndices, keep, combine_filters);
            } else {
                passingIndices = keep;
            }
        }
    }

    return passingIndices;
}

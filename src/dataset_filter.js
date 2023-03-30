import {isArray} from 'lodash';

import {getVarNameType} from './VectorUtil';

function combine(a, b, op) {
  return op === 'or'
    ? new Set([...a, ...b])
    : new Set([...a].filter((x) => b.has(x)));
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

function createSingleFilterFunction(filter) {
  let {operation, value} = filter;
  let applyFunction;
  if (operation === 'in') {
    value = new Set(value);
    applyFunction = (d) => value.has(d);
  } else if (operation === '>') {
    applyFunction = (d) => d > value;
  } else if (operation === '=') {
    applyFunction = (d) => d === value;
  } else if (operation === '<') {
    applyFunction = (d) => d < value;
  } else if (operation === '!=') {
    applyFunction = (d) => d !== value;
  } else if (operation === '>=') {
    applyFunction = (d) => d >= value;
  } else if (operation === '<=') {
    applyFunction = (d) => d <= value;
  } else {
    throw new Error('Unknown filter: ' + operation);
  }
  return applyFunction;
}

export function createFilterFunction(filter) {
  let {operation, value, invert} = filter;
  let applyFunction;
  if (isArray(operation)) {
    const functions = [];
    for (let i = 0; i < operation.length; i++) {
      const f = createSingleFilterFunction({
        operation: operation[i],
        value: value[i],
      });
      functions.push(f);
    }
    applyFunction = (d) => {
      for (let i = 0; i < functions.length; i++) {
        if (!functions[i](d)) {
          return invert;
        }
      }
      return !invert;
    };
  } else {
    applyFunction = createSingleFilterFunction(filter);
  }
  return applyFunction;
}

export function getPassingFilterIndices(cachedData, datasetFilter) {
  let passingIndices = null;
  if (datasetFilter) {
    let userFilters = datasetFilter.filters || [];
    let combineFilters = datasetFilter.combine || 'and';
    for (let i = 0; i < userFilters.length; i++) {
      const filterObject = userFilters[i];
      const {field, value} = filterObject;
      let keep = null;

      if (field === '__index') {
        // ['__index', 'in', Array(1218)] for lasso or box filters
        keep = new Set(value);
      } else {
        const nameType = getVarNameType(field);
        const series = cachedData[nameType.name];
        keep = getIndices(series, createFilterFunction(filterObject));
      }
      if (passingIndices) {
        passingIndices = combine(passingIndices, keep, combineFilters);
      } else {
        passingIndices = keep;
      }
    }
  }
  return passingIndices;
}

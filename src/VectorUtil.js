import natsort from 'natsort';
import {SlicedVector} from './SlicedVector';


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

export function getIndices(v, f) {
    let result = new Set();
    for (let i = 0, size = v.size(); i < size; i++) {
        if (f(v.get(i))) {
            result.add(i);
        }
    }
    return result;
}

export function groupedStats(dimensions, varMeasures) {
    let allResults = [];
    dimensions.forEach(v => {

        let valueToIndices = {};
        for (let i = 0, size = v.size(); i < size; i++) {
            let value = v.get(i);
            let indices = valueToIndices[value];
            if (indices === undefined) {
                indices = [];
                valueToIndices[value] = indices;
            }
            indices.push(i);
        }
        let categories = Object.keys(valueToIndices);
        categories.sort(natsort());
        const results = {'categories': categories, 'name': v.getName(), 'values': []};
        varMeasures.forEach(v => {
            results.values.push({fractionExpressed: [], mean: [], name: v.getName()});
        });
        categories.forEach(category => {
            let indices = valueToIndices[category];
            varMeasures.forEach((v, index) => {
                let slicedVector = new SlicedVector(v, indices);
                let xStats = stats(slicedVector);
                results.values[index].mean.push(xStats.mean);
                results.values[index].fractionExpressed.push(xStats.numExpressed / slicedVector.size());
            });
        });

        allResults.push(results);
    });
    return allResults;
}

export function stats(v) {
    let min = Number.MAX_VALUE;
    let max = -Number.MAX_VALUE;
    let sum = 0;
    let numExpressed = 0;

    for (let i = 0, size = v.size(); i < size; i++) {
        const value = v.get(i);
        min = value < min ? value : min;
        max = value > max ? value : max;
        sum += value;
        if (value !== 0) {
            numExpressed++;
        }
    }
    return {
        'min': min, 'max': max, 'sum': sum, 'mean': sum / v.size(), 'numExpressed': numExpressed,
    };
}




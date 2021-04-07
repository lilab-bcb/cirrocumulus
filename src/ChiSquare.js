import {regularizedGammaPValue} from './RegularizedGamma';

function cumulativeProb(x, shape, scale) {
    if (x <= 0) {
        return 0;
    } else if (x >= Number.POSITIVE_INFINITY) {
        return 1;
    }

    return regularizedGammaPValue(shape, x / scale);
}

export function chiSquare(counts) {
    const rowSums = [];

    const nrows = counts.length;
    const ncols = counts[0].length;
    let total = 0;
    for (let i = 0; i < nrows; i++) {
        let sum = 0;
        for (let j = 0; j < ncols; j++) {
            sum += counts[i][j];
        }
        rowSums.push(sum);
        total += sum;
    }
    const columnSums = [];
    for (let j = 0; j < ncols; j++) {
        let sum = 0;
        for (let i = 0; i < nrows; i++) {
            sum += counts[i][j];
        }
        columnSums.push(sum);
    }
    const expectedTable = []; //   (row sum * column sum)/N
    for (let i = 0; i < nrows; i++) {
        const expected = [];
        for (let j = 0; j < ncols; j++) {
            expected.push(rowSums[i] * columnSums[j] / total);
        }
        expectedTable.push(expected);
    }
    let statistic = 0;
    for (let i = 0; i < nrows; i++) {
        for (let j = 0; j < ncols; j++) {
            statistic += Math.pow(counts[i][j] - expectedTable[i][j], 2) / expectedTable[i][j];
        }
    }
    const degreesOfFreedom = (nrows - 1) * (ncols - 1);
    const p = 1 - cumulativeProb(statistic, degreesOfFreedom / 2, 2);
    return {expected: expectedTable, statistic: statistic, p: p, df: degreesOfFreedom};
}




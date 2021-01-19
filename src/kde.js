import {quantileSorted} from 'd3-array';


// Scott, D. W. (1992) Multivariate Density Estimation: Theory, Practice, and
// Visualization. Wiley.
// export function nrd(x) {
//     let h = iqr(x) / 1.34;
//     return 1.06 * Math.min(Math.sqrt(variance(x)), h)
//         * Math.pow(x.length, -1 / 5);
// }

// Silverman, B. W. (1986) Density Estimation. London: Chapman and Hall.
// export function nrd0(x) {
//     let hi = Math.sqrt(variance(x));
//     let lo;
//     if (!(lo = Math.min(hi, iqr(x) / 1.34))) {
//         (lo = hi) || (lo = Math.abs(x[1])) || (lo = 1);
//     }
//     return .9 * lo * Math.pow(x.length, -.2);
// }

export function nrd0(stats) {
    let hi = Math.sqrt(stats.variance);
    const iqr = stats.q3 - stats.q1;
    let lo;
    if (!(lo = Math.min(hi, iqr / 1.34))) {
        (lo = hi) || (lo = Math.abs(stats.x0)) || (lo = 1);
    }
    return .9 * lo * Math.pow(stats.n, -.2);
}

export function boxplotStats(x) {
    x = x.slice().sort((a, b) => a - b);
    const q3 = quantileSorted(x, 0.75);
    const q1 = quantileSorted(x, 0.25);
    const q50 = quantileSorted(x, 0.5);
    const w = 1.5;
    let upperAdjacentValue = -Number.MAX_VALUE;
    let lowerAdjacentValue = Number.MAX_VALUE;
    // The upper adjacent value (UAV) is the largest observation that is
    // less than or equal to
    // the upper inner fence (UIF), which is the third quartile plus
    // 1.5*IQR.
    //
    // The lower adjacent value (LAV) is the smallest observation that is
    // greater than or equal
    // to the lower inner fence (LIF), which is the first quartile minus
    // 1.5*IQR.
    let upperOutlier = q3 + w * (q3 - q1);
    let lowerOutlier = q1 - w * (q3 - q1);
    let sum = 0;
    for (let i = 0, n = x.length; i < n; i++) {
        const value = x[i];
        if (value <= upperOutlier) {
            upperAdjacentValue = Math.max(upperAdjacentValue, value);
        }
        if (value >= lowerOutlier) {
            lowerAdjacentValue = Math.min(lowerAdjacentValue, value);
        }
        sum += value;
    }
    const mean = sum / x.length;
    if (lowerAdjacentValue > q1) {
        lowerAdjacentValue = q1;
    }
    if (upperAdjacentValue < q3) {
        upperAdjacentValue = q3;
    }
    return {
        mean: mean,
        upperAdjacentValue: upperAdjacentValue,
        lowerAdjacentValue: lowerAdjacentValue,
        median: q50,
        q3: q3,
        q1: q1,
        variance: variance(x, mean),
        n: x.length,
        x0: x[0]
    };
}

function variance(x, mean) {
    let sum = 0;
    let n = x.length;

    for (let i = 0; i < n; i++) {
        const v = x[i];
        let diff = v - mean;
        diff = diff * diff;
        sum += diff;
    }
    if (n <= 1) {
        return NaN;
    }
    n = n - 1;
    if (n < 1) {
        n = 1;
    }
    return sum / n;
}


function gaussian(v) {
    return (1 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * v * v);
}

function createKDE(bandwidth, kernel, values) {
    const len = values.length;
    const factor = 1 / (len * bandwidth);

    return function (x) {
        let sum = 0;
        for (let i = 0; i < len; i++) {
            sum += kernel((x - values[i]) / bandwidth);
        }
        return factor * sum;
    };
}


export function density(values, bandwidth, gridsize = 200) {
    let min = Number.MAX_VALUE;
    let max = -Number.MAX_VALUE;
    for (let i = 0, n = values.length; i < n; i++) {
        const value = values[i];
        min = Math.min(min, value);
        max = Math.max(max, value);
    }
    const span = max - min;
    const step = span / (gridsize - 1);
    const x = new Float32Array(gridsize);
    const y = new Float32Array(gridsize);
    let maxKDE = 0;
    const kde = createKDE(bandwidth, gaussian, values);
    for (let i = 0, j = min; i < gridsize; i++, j += step) {
        const v = kde(j);
        x[i] = j;
        y[i] = v;
        maxKDE = Math.max(maxKDE, v);
    }
    return {x: x, y: y, max: maxKDE};
}


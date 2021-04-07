import {regularizedGammaQValue} from './RegularizedGamma';
import {rankdata} from './util';

const SQRT2 = Math.sqrt(2.0);
const EXTREME_VALUE_BOUND = 40;

function calculateAsymptoticPValue(Umin, n1, n2) {

    /* long multiplication to avoid overflow (double not used due to efficiency
     * and to avoid precision loss)
     */
    const n1n2prod = n1 * n2;

    // http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U#Normal_approximation
    const EU = n1n2prod / 2.0;
    const VarU = n1n2prod * (n1 + n2 + 1) / 12.0;

    const z = (Umin - EU) / Math.sqrt(VarU);

    return 2 * standardNormalCumulativeProbability(z, 0, 1);
}

function erfc(x) {
    if (Math.abs(x) > EXTREME_VALUE_BOUND) {
        return x > 0 ? 0 : 2;
    }
    const ret = regularizedGammaQValue(0.5, x * x, 1e-15, 10000);
    return x < 0 ?
        2 - ret :
        ret;
}

function standardNormalCumulativeProbability(x, mean, standardDeviation) {
    const dev = x - mean;
    if (Math.abs(dev) > 40 * standardDeviation) {
        return dev < 0 ? 0.0 : 1.0;
    }
    return 0.5 * erfc(-dev / (standardDeviation * SQRT2));

}

export function mannWhitney(x, y) {
    const Umax = mannWhitneyU(x, y);

    /*
     * It can be shown that U1 + U2 = n1 * n2
     */
    const Umin = x.length * y.length - Umax;
    const p = calculateAsymptoticPValue(Umin, x.length, y.length);
    return {statistic: Umax, p: p};
}

function mannWhitneyU(x, y) {
    const z = []; // combine x and y
    for (let i = 0, n = x.length; i < n; i++) {
        z.push(x[i]);
    }
    for (let i = 0, n = y.length; i < n; i++) {
        z.push(y[i]);
    }
    const ranks = rankdata(z);
    let sumRankX = 0;

    /*
     * The ranks for x is in the first x.length entries in ranks because x
     * is in the first x.length entries in z
     */
    for (let i = 0; i < x.length; ++i) {
        sumRankX += ranks[i];
    }

    /*
     * U1 = R1 - (n1 * (n1 + 1)) / 2 where R1 is sum of ranks for sample 1,
     * e.g. x, n1 is the number of observations in sample 1.
     */
    const U1 = sumRankX - (x.length * (x.length + 1)) / 2;

    /*
     * It can be shown that U1 + U2 = n1 * n2
     */
    const U2 = x.length * y.length - U1;

    return Math.max(U1, U2);
}

import {color} from 'd3-color';
import {scaleLinear} from 'd3-scale';
import * as scaleChromatic from 'd3-scale-chromatic';
import simplify from 'simplify-js';
import {getColors} from './ThreeUtil';

export const interpolators = {};
interpolators['Diverging'] = [
    'interpolateBrBG',
    'interpolatePRGn',
    'interpolatePiYG',
    'interpolatePuOr',
    'interpolateRdBu',
    'interpolateRdGy',
    'interpolateRdYlBu',
    'interpolateRdYlGn',
    'interpolateSpectral'];

interpolators['Sequential (Single Hue)'] = [
    'interpolateBlues',
    'interpolateGreens',
    'interpolateGreys',
    'interpolateOranges',
    'interpolatePurples',
    'interpolateReds'];

interpolators['Sequential (Multi-Hue)'] = [
    'interpolateViridis',
    'interpolateInferno',
    'interpolateMagma',
    'interpolatePlasma',
    'interpolateWarm',
    'interpolateCool',
    'interpolateCubehelixDefault',
    'interpolateBuGn',
    'interpolateBuPu',
    'interpolateGnBu',
    'interpolateOrRd',
    'interpolatePuBuGn',
    'interpolatePuBu',
    'interpolatePuRd',
    'interpolateRdPu',
    'interpolateYlGnBu',
    'interpolateYlGn',
    'interpolateYlOrBr',
    'interpolateYlOrRd'];

interpolators['Cyclical'] = ['interpolateRainbow', 'interpolateSinebow'];


// const TWENTY_COLORS = [
//     '#1f77b4', '#aec7e8', '#ff7f0e',
//     '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd',
//     '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
//     '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5'];

export const CATEGORY_20B = [
    '#393b79', '#5254a3', '#6b6ecf',
    '#9c9ede', '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31',
    '#bd9e39', '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b',
    '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6'];
export const CATEGORY_20C = [
    '#3182bd', '#6baed6', '#9ecae1',
    '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2', '#31a354',
    '#74c476', '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc',
    '#dadaeb', '#636363', '#969696', '#bdbdbd', '#d9d9d9'];

export function getChartSize() {
    // leave room for drawer and header
    return {width: window.screen.availWidth - 280, height: window.screen.availHeight - 220};
}

/**
 *
 * @param array. Array of format,data
 */
export function setClipboardData(clipboardData) {
    const container = document.activeElement;
    const isRTL = document.documentElement.getAttribute('dir') == 'rtl';
    const fakeElem = document.createElement('div');
    fakeElem.contentEditable = true;

    // Prevent zooming on iOS
    fakeElem.style.fontSize = '12pt';
    // Reset box model

    fakeElem.style.border = '0';
    fakeElem.style.padding = '0';
    fakeElem.style.margin = '0';
    // Move element out of screen horizontally
    fakeElem.style.position = 'absolute';
    fakeElem.style[isRTL ? 'right' : 'left'] = '-999999px';
    // Move element to the same position vertically
    fakeElem.style.top = (window.pageYOffset || document.documentElement.scrollTop) + 'px';
    fakeElem.setAttribute('readonly', '');
    //fakeElem.innerHTML = html;
    const copyListener = (e) => {

        clipboardData.forEach(function (elem) {
            e.clipboardData.setData(elem.format, elem.data);
        });

        e.preventDefault();
        e.stopPropagation();
        e.stopImmediatePropagation();
        fakeElem.removeEventListener('copy', copyListener);
    };
    fakeElem.addEventListener('copy', copyListener);

    container.appendChild(fakeElem);

    const selection = window.getSelection();
    const range = document.createRange();
    range.selectNodeContents(fakeElem);
    selection.removeAllRanges();
    selection.addRange(range);
    const fakeHandlerCallback = (event) => {
        container.removeChild(fakeElem);
        container.removeEventListener('click', fakeHandlerCallback);
    };
    document.execCommand('copy');
    container.addEventListener('click', fakeHandlerCallback);
};

export function updateTraceColors(traceInfo) {
    if (traceInfo.isImage) {
        let colors = [];
        let colorScale = traceInfo.colorScale;
        const colorMapper = rgb => rgb.formatHex();
        for (let i = 0, n = traceInfo.npoints; i < n; i++) {
            let rgb = color(colorScale(traceInfo.values[i]));
            colors.push(colorMapper(rgb));
        }
        traceInfo.colors = colors;
    } else {
        traceInfo.colors = getColors(traceInfo);
    }
}


/**
 * Computes the rank using the given index array. The index array can be
 * obtained from the indexSort method. Does not handle ties.
 *
 */
function rankIndexArray(index) {
    const rank = [];
    const n = index.length;
    for (let j = 0; j < n; j++) {
        rank[index[j]] = j + 1;
    }
    return rank;
};

function indexSort(array, ascending) {
    const pairs = [];
    for (let i = 0, length = array.length; i < length; i++) {
        pairs.push({
            value: array[i],
            index: i
        });
    }
    return indexSortPairs(pairs, ascending);
};

function indexSortPairs(array, ascending) {
    if (ascending) {
        array.sort(function (a, b) {
            return (a.value < b.value ? -1 : (a.value === b.value ? (a.index < b.index ? -1 : 1) : 1));
        });
    } else {
        array.sort(function (a, b) {
            return (a.value < b.value ? 1 : (a.value === b.value ? (a.index < b.index ? 1 : -1) : -1));
        });
    }
    const indices = [];
    array.forEach(function (item) {
        indices.push(item.index);
    });
    return indices;
};

/**
 * Computes the False Discovery Rate using the BH procedure.
 *
 * @param nominalPValues
 *            Array of nominal p-values.
 */
export function fdr(nominalPValues) {
    const size = nominalPValues.length;
    const pValueIndices = indexSort(nominalPValues, true);
    const ranks = rankIndexArray(pValueIndices);

    // check for ties
    for (let i = pValueIndices.length - 1; i > 0; i--) {
        const bigPValue = nominalPValues[pValueIndices[i]];
        const smallPValue = nominalPValues[pValueIndices[i - 1]];
        if (bigPValue === smallPValue) {
            ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
        }
    }

    const fdr = new Float32Array(size);
    for (let i = 0; i < size; i++) {
        const rank = ranks[i];
        const p = nominalPValues[i];
        fdr[i] = (p * size) / rank;
    }

    // ensure fdr is monotonically decreasing
    const pIndices = indexSort(nominalPValues, false);
    for (let i = 0; i < pIndices.length - 1; i++) {
        const highIndex = pIndices[i];
        const lowIndex = pIndices[i + 1];
        fdr[lowIndex] = Math.min(fdr[lowIndex], fdr[highIndex]);
    }
    for (let i = 0; i < size; i++) {
        fdr[i] = Math.min(fdr[i], 1);
    }
    return fdr;
};

export function isPointInside(point, vs) {
    // ray-casting algorithm based on
    // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    const x = point.x, y = point.y;

    let inside = false;
    for (let i = 0, j = vs.length - 1; i < vs.length; j = i++) {
        const xi = vs[i].x, yi = vs[i].y;
        const xj = vs[j].x, yj = vs[j].y;

        const intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }

    return inside;
}


export function arrayToSvgPath(lassoPathArray) {
    if (lassoPathArray.length > 1) {
        lassoPathArray = simplify(lassoPathArray);
    }
    let svgPath = 'M ' + lassoPathArray[0].x + ' ' + lassoPathArray[0].y;
    for (let i = 1; i < lassoPathArray.length; i++) {
        svgPath += ' L ' + lassoPathArray[i].x + ' ' + lassoPathArray[i].y;
    }
    svgPath += ' Z';
    return svgPath;
}

export function getRgbScale() {
    return scaleLinear().domain([0, 255]).range([0, 1]);
}

export function fixInterpolatorName(name) {
    if (!name.startsWith("interpolate")) {
        name = "interpolate" + name;
    }
    return name;
}

export function getInterpolator(name) {
    return scaleChromatic[fixInterpolatorName(name)];
}

export function convertPointsToBins(points, allBins) {
    let bins = [];
    for (let i = 0, n = points.length; i < n; i++) {
        bins.push(allBins[points[i]]);
    }
    return bins;
}

export function convertBinsToPoints(bins, selectedBins) {
    let points = [];
    let binIndex = 0;

    for (let i = 0, n = bins.length; i < n; i++) {
        if (bins[i] === selectedBins[binIndex]) {
            points.push(i);
            binIndex++;
        }
    }

    // selectedBins = new Set(selectedBins);
    // for (let i = 0, n = bins.length; i < n; i++) {
    //     if (selectedBins.has(bins[i])) {
    //         points.push(i);
    //     }
    // }

    return points;
}


export function splitSearchTokens(tokens) {
    let X = [];
    let obs = [];
    let obsCat = [];
    let featureSets = [];
    tokens.forEach(token => {
        if (token.type === 'X') {
            X.push(token.value);
        } else if (token.type === 'obs') {
            obs.push(token.value);
        } else if (token.type === 'obsCat') {
            obsCat.push(token.value);
        } else if (token.type === 'featureSet') {
            featureSets.push(token.value);
        } else {
            console.log('Unknown type: ' + token);
        }
    });
    return {X: X, obs: obs, obsCat: obsCat, featureSets: featureSets};
}
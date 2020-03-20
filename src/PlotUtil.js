import {scaleLinear} from 'd3-scale';
import * as scaleChromatic from 'd3-scale-chromatic';
import simplify from 'simplify-js';

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
    return {width: window.screen.availWidth - 300, height: window.screen.availHeight - 200};
}

export function isPointInside(point, vs) {
    // ray-casting algorithm based on
    // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    var x = point.x, y = point.y;

    var inside = false;
    for (var i = 0, j = vs.length - 1; i < vs.length; j = i++) {
        var xi = vs[i].x, yi = vs[i].y;
        var xj = vs[j].x, yj = vs[j].y;

        var intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }

    return inside;
}

export function drawScatter2d(context, chartSize, traceInfo, markerSize, markerOpacity, unselectedMarkerOpacity, selection, color) {
    let height = chartSize.height;
    let width = chartSize.width;

    let xmin = Number.MAX_VALUE;
    let xmax = -Number.MAX_VALUE;
    let ymin = Number.MAX_VALUE;
    let ymax = -Number.MAX_VALUE;
    for (let i = 0, n = traceInfo.x.length; i < n; i++) {
        let x = traceInfo.x[i];
        let y = traceInfo.y[i];
        xmin = x < xmin ? x : xmin;
        xmax = x > xmax ? x : xmax;
        ymin = y < ymin ? y : ymin;
        ymax = y > ymax ? y : ymax;
    }
    const xToPixScale = scaleLinear().domain([xmin, xmax]).range([markerSize, width - markerSize]);
    const yToPixScale = scaleLinear().domain([ymin, ymax]).range([height - markerSize, markerSize]);
    const PI2 = 2 * Math.PI;
    const colorScale = scaleLinear().domain([0, 1]).range([0, 255]);
    for (let i = 0, n = traceInfo.x.length; i < n; i++) {
        const isSelected = selection.size === 0 || selection.has(i);
        const x = traceInfo.x[i];
        const y = traceInfo.y[i];
        const xpix = xToPixScale(x);
        const ypix = yToPixScale(y);
        const c = color[i];
        const alpha = isSelected ? markerOpacity : unselectedMarkerOpacity;
        context.fillStyle = 'rgba(' + colorScale(c.r) + ',' + colorScale(c.g) + ',' + colorScale(c.b) + ',' + alpha + ')';
        context.beginPath();
        context.arc(xpix, ypix, markerSize, 0, PI2);
        context.closePath();
        context.fill();
    }
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

class PlotUtil {


    static convertPointsToBins(points, allBins) {
        let bins = [];
        for (let i = 0, n = points.length; i < n; i++) {
            bins.push(allBins[points[i]]);
        }
        return bins;
    }

    static convertBinsToPoints(bins, selectedBins) {
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

    static createPlotConfig() {
        return {
            showLink: false,
            responsive: false,
            displaylogo: false,
            scrollZoom: false,
            displayModeBar: true,
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines', 'sendDataToCloud']
        };
    }

    static createDotPlotConfig() {
        return {
            showLink: false,
            scrollZoom: false,
            responsive: false,
            displaylogo: false,
            displayModeBar: true,
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines', 'sendDataToCloud', 'zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d']
        };
    }


    static createDotPlotAxis() {
        return {
            showbackground: false,
            autorange: false,
            fixedrange: true,
            showgrid: true,
            zeroline: false,
            showline: false,
            title: '',
            // type: 'category'
            autotick: false
        };
    }

    static createDotPlotLayout(options) {
        let {width, height} = options;
        let layout = {
            hovermode: 'closest',
            dragmode: 'select',
            width: width,
            height: height,
            margin: {
                l: 0,
                b: 0,
                r: 0,
                t: 0,
                autoexpand: false
            },
            fixedrange: true,
            legend: {yanchor: 'top'},
            autosize: true,
            displaylogo: false,
            showlegend: false,
            font: {
                family: 'Roboto Condensed,Helvetica,Arial,sans-serif'
            }
        };

        layout.xaxis = PlotUtil.createDotPlotAxis();
        layout.yaxis = PlotUtil.createDotPlotAxis();
        return layout;
    }
}

export default PlotUtil;

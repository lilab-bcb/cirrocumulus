import * as scaleChromatic from 'd3-scale-chromatic';

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

export function isPlotlyBug(el, newTrace) {
    // https://github.com/plotly/plotly.js/issues/3405
    const threshold = 100000;
    const oldData = el.data;
    const oldSize = oldData != null ? oldData[0].x.length : 0;
    const newSize = newTrace.x != null ? newTrace.x.length : 0;
    return ((oldSize > threshold) && (newSize <= threshold));
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
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines', 'sendDataToCloud']
        };
    }

    static createDotPlotConfig() {
        return {
            showLink: false,
            responsive: false,
            displaylogo: false,
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines', 'sendDataToCloud', 'zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d']
        };
    }


    static createEmbeddingAxis() {
        return {
            showbackground: false,
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            title: '',
            autotick: false,
            ticks: '',
            showticklabels: false,

        };
    }


    static getEmbeddingChartSize(size) {
        let maxSize = Math.floor(Math.min(window.screen.availWidth - 240, window.screen.availHeight - 190));
        return Math.floor(maxSize / size);
    }

    static createEmbeddingLayout(options) {
        let {size, is3d} = options;
        size = PlotUtil.getEmbeddingChartSize(size);

        let layout = {
            hovermode: 'closest',
            dragmode: 'select',
            width: size,
            height: size,
            margin: {
                l: 0,
                b: 0,
                r: 0,
                t: 0,
                autoexpand: false
            },
            legend: {yanchor: 'top'},
            autosize: false,
            displaylogo: false,
            showlegend: false,

        };
        if (is3d) {
            layout.scene = {
                xaxis: PlotUtil.createEmbeddingAxis(),
                yaxis: PlotUtil.createEmbeddingAxis(),
                zaxis: PlotUtil.createEmbeddingAxis(),
            };
        } else {
            layout.xaxis = PlotUtil.createEmbeddingAxis();
            layout.yaxis = PlotUtil.createEmbeddingAxis();
        }
        return layout;
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

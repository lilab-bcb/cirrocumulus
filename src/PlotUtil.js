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

    static createDotPlotAxis() {
        return {
            showbackground: false,
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            title: '',
            type: 'category'
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
                autoexpand: true
            },
            legend: {yanchor: 'top'},
            autosize: true,
            displaylogo: false,
            showlegend: false,
        };

        layout.xaxis = PlotUtil.createDotPlotAxis();
        layout.yaxis = PlotUtil.createDotPlotAxis();
        return layout;
    }
}

export default PlotUtil;

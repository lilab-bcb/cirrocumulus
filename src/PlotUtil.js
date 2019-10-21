import {format} from 'd3-format';
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
const intFormat = format(',');
const percentFormat = format('.1f');


export function getLegendSizeHelper(selectedCountMap, scale, index, selectionCount) {
    if (scale.valueCounts.total == null) {
        // set total lazily
        let total = 0;
        for (let i = 0, n = scale.valueCounts.counts.length; i < n; i++) {
            total += scale.valueCounts.counts[i];
        }
        scale.valueCounts.total = total;
    }

    let count = scale.valueCounts.counts[index];
    let total = scale.valueCounts.total;
    let percent = count / total;
    let percentSelected = Number.NaN;
    let title = intFormat(count) + ' / ' + intFormat(total) + ' (' + percentFormat(100 * percent) + '% of total)';
    let selectionTitle;
    if (selectedCountMap != null) {
        let d = scale.valueCounts.values[index];
        let selectedCount = selectedCountMap[d] || 0;
        percentSelected = selectedCount / selectionCount;
        selectionTitle = intFormat(selectedCount) + ' / ' + intFormat(selectionCount) + ' % of selection (' + percentFormat(100 * percentSelected) + '%)';
    }

    return {
        percentTotal: percent,
        percentSelected: percentSelected,
        title: title,
        total: intFormat(total),
        selectionTitle: selectionTitle
    };
}


export function getInterpolator(name) {
    if (!name.startsWith("interpolate")) {
        name = "interpolate" + name;
    }
    return scaleChromatic[name];
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
            showLink: true,
            responsive: false,
            displaylogo: false,
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines'],// 'sendDataToCloud'
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
                l: 50,
                b: 50,
                r: 100,
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

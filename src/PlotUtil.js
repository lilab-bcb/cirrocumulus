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

export function getInterpolator(name) {
    if (!name.startsWith("interpolate")) {
        name = "interpolate" + name;
    }
    return scaleChromatic[name];
}

class PlotUtil {

    static createPlotConfig() {
        return {
            showLink: true,
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines'],// 'sendDataToCloud'
        };
    }

    static fixLegend(node) {
        requestAnimationFrame(() => {
            if (node) {
                const scatterpts = node.querySelectorAll('.scatterpts');
                if (scatterpts) {
                    scatterpts.forEach(pt => {
                        pt.setAttribute('d', 'M6,0A6,6 0 1,1 0,-6A6,6 0 0,1 6,0Z');
                        pt.style.opacity = 1;
                    });
                }
            }
        });
    }

    static createAxis(options) {
        return options.embedding ? {
            showbackground: false,
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            autotick: false,
            ticks: '',
            showticklabels: false,
            title: '',
        } : {
            showbackground: false,
            autorange: true,
            showgrid: false,
            zeroline: false,
            showline: false,
            title: '',
        };
    }

    static createPlotLayout(options) {
        let legendWidth = options.legend || 0;
        let size = Math.floor(Math.min(window.screen.availWidth - 220, window.screen.availHeight * .95));
        let layout = {
            hovermode: 'closest',
            dragmode: 'select',
            title: {text: options.title, font: {size: 12}},
            width: size + legendWidth,
            height: size,
            margin: {
                l: !options.embedding ? 50 : 0,
                b: !options.embedding ? 50 : 0,
                r: legendWidth + (!options.embedding ? 100 : 0),
                t: options.title ? 22 : 0,
                autoexpand: true
            },
            legend: {yanchor: 'top'},
            autosize: !options.embedding,
            displaylogo: false,
            showlegend: options.legend > 0 ? true : false,
        };
        if (options.is3d) {
            layout.scene = {
                xaxis: PlotUtil.createAxis(options),
                yaxis: PlotUtil.createAxis(options),
                zaxis: PlotUtil.createAxis(options),
            };
        } else {
            layout.xaxis = PlotUtil.createAxis(options);
            layout.yaxis = PlotUtil.createAxis(options);
        }
        return layout;
    }
}

export default PlotUtil;

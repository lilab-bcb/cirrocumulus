import {color} from 'd3-color';
import {makeStyles, ScatterPlot, ScatterPlotVisualizerSprites} from 'scatter-gl';
import {getRgbScale} from './PlotUtil';

function scaleLinear(value, domain, range) {
    const domainDifference = domain[1] - domain[0];
    const rangeDifference = range[1] - range[0];

    const percentDomain = (value - domain[0]) / domainDifference;
    return percentDomain * rangeDifference + range[0];
}

export function createScatterPlot(containerElement) {
    const styles = makeStyles();

    const scatterPlot = new ScatterPlot(containerElement, {
        camera: {},
        selectEnabled: false,
        styles: styles,

    });

    const activeVisualizers = [];
    const visualizer = new ScatterPlotVisualizerSprites(styles);
    activeVisualizers.push(visualizer);
    scatterPlot.setActiveVisualizers(activeVisualizers);
    scatterPlot.setInteractionMode('SELECT');
    return scatterPlot;
}

export function getColors(trace) {
    const RGBA_NUM_ELEMENTS = 4;
    const rgbScale = getRgbScale();

    let dst = 0;
    let colorScale = trace.colorScale;
    const n = trace.x.length;
    const colors = new Float32Array(n * RGBA_NUM_ELEMENTS);

    for (let i = 0; i < n; ++i) {
        let c = color(colorScale(trace.values[i]));
        colors[dst++] = rgbScale(c.r);
        colors[dst++] = rgbScale(c.g);
        colors[dst++] = rgbScale(c.b);
        colors[dst++] = 1;
    }
    return colors;
}

export function getPositions(trace) {
    const SCATTER_PLOT_CUBE_LENGTH = 2;
    let xExtent = [Infinity, -Infinity];
    let yExtent = [Infinity, -Infinity];
    let zExtent = [Infinity, -Infinity];
    const npoints = trace.x.length;
    const is3d = trace.z != null;
    // Determine max and min of each axis of our data.
    for (let i = 0; i < npoints; i++) {
        const x = trace.x[i];
        if (x < xExtent[0]) xExtent[0] = x;
        if (x > xExtent[1]) xExtent[1] = x;

        const y = trace.y[i];
        if (y < yExtent[0]) yExtent[0] = y;
        if (y > yExtent[1]) yExtent[1] = y;
        if (is3d) {
            const z = trace.z[i];
            if (z < zExtent[0]) zExtent[0] = z;
            if (z > zExtent[1]) zExtent[1] = z;
        }
    }

    const getRange = (extent) => Math.abs(extent[1] - extent[0]);
    const xRange = getRange(xExtent);
    const yRange = getRange(yExtent);
    const zRange = is3d ? getRange(zExtent) : 0;
    const maxRange = Math.max(xRange, yRange);
    const halfCube = SCATTER_PLOT_CUBE_LENGTH / 2;
    const makeScaleRange = (range, base) => [
        -base * (range / maxRange),
        base * (range / maxRange),
    ];
    const xScale = makeScaleRange(xRange, halfCube);
    const yScale = makeScaleRange(yRange, halfCube);
    const zScale = makeScaleRange(zRange, halfCube);

    const positions = new Float32Array(npoints * 3);
    let dst = 0;

    for (let i = 0; i < npoints; i++) {

        positions[dst++] = scaleLinear(trace.x[i], xExtent, xScale);
        positions[dst++] = scaleLinear(trace.y[i], yExtent, yScale);

        if (is3d) {
            positions[dst++] = scaleLinear(trace.z[i], zExtent, zScale);
        } else {
            positions[dst++] = 0.0;
        }
    }
    return positions;

}

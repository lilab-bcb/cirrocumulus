import {color} from 'd3-color';
import {makeStyles, ScatterPlot, ScatterPlotVisualizerCanvasLabels, ScatterPlotVisualizerSprites} from 'scatter-gl';
import {Color} from 'three';
import {getVisualizer, LABELS_VISUALIZER_ID, POINT_VISUALIZER_ID} from './ScatterChartThree';
import {getRgbScale} from './util';

function scaleLinear(value, domain, range) {
    const domainDifference = domain[1] - domain[0];
    const rangeDifference = range[1] - range[0];

    const percentDomain = (value - domain[0]) / domainDifference;
    return percentDomain * rangeDifference + range[0];
}


export function createScatterPlot(containerElement, premultipliedAlpha, labels) {
    const styles = makeStyles();
    styles.label3D.fontSize = 40;

    const scatterPlot = new ScatterPlot(containerElement, {
        camera: {},
        selectEnabled: false,
        styles: styles,

    }, premultipliedAlpha); // toDataUrl images are flipped on Safari when premultipliedAlpha is false
    let visualizers = [new ScatterPlotVisualizerSprites(styles)];
    if (labels) {
        visualizers.push(new ScatterPlotVisualizerCanvasLabels(containerElement, styles));
    }
    scatterPlot.setActiveVisualizers(visualizers);
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
    const zRange = is3d ? getRange(zExtent) : 1;
    const maxRange = Math.max(xRange, yRange, zRange);
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

export function getCategoryLabelsPositions(traceInfo, categoricalNames) {
    const categoryToPosition = {};
    let ncategories = 0;
    const isImage = traceInfo.isImage;
    for (let i = 0, j = 0; i < traceInfo.npoints; i++, j += 3) {
        let value = traceInfo.values[i];
        let p = categoryToPosition[value];
        if (p === undefined) {
            p = {count: 0, position: [0, 0, 0]};
            categoryToPosition[value] = p;
            ncategories++;
        }
        p.count++;
        if (isImage) {
            p.position[0] += traceInfo.x[i];
            p.position[1] += traceInfo.y[i];
        } else {
            p.position[0] += traceInfo.positions[j];
            p.position[1] += traceInfo.positions[j + 1];
            p.position[2] += traceInfo.positions[j + 2];
        }
    }
    let labelStrings = [];
    let labelPositions = new Float32Array(ncategories * 3);
    let positionIndex = 0;
    let categoryObject = categoricalNames[traceInfo.name];
    if (categoryObject === undefined) {
        categoryObject = {};
    }
    for (let category in categoryToPosition) {
        let renamedCategory = categoryObject[category];
        if (renamedCategory !== undefined) {
            labelStrings.push(renamedCategory);
        } else {
            labelStrings.push(category);
        }
        let p = categoryToPosition[category];
        labelPositions[positionIndex] = p.position[0] / p.count;
        labelPositions[positionIndex + 1] = p.position[1] / p.count;
        labelPositions[positionIndex + 2] = p.position[2] / p.count;
        positionIndex += 3;
    }

    return {labels: labelStrings, positions: labelPositions};
}

export function updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize, showLabels = false, categoricalNames = {}, showFog = true, showAxis = true, darkMode = false) {
    const colors = traceInfo.colors;
    const positions = traceInfo.positions;
    const is3d = traceInfo.z != null;
    for (let i = 0, j = 3, k = 2; i < traceInfo.npoints; i++, j += 4, k += 3) {
        const isSelected = selection.size === 0 || selection.has(i);
        colors[j] = isSelected ? markerOpacity : unselectedMarkerOpacity;
        if (!is3d) {
            positions[k] = isSelected ? 1 : 0;
        }
    }
    scatterPlot.scene.background = darkMode ? new Color("rgb(0, 0, 0)") : null;
    scatterPlot.setDimensions(traceInfo.dimensions);
    let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
    spriteVisualizer.styles.fog.enabled = showFog;
    const axes = scatterPlot.scene.getObjectByName('axes');
    if (axes) {
        axes.visible = showAxis;
    }
    scatterPlot.setPointColors(colors);
    scatterPlot.setPointPositions(positions);

    // const {scaleDefault, scaleSelected, scaleHover} = scatterPlot.styles.point;

    const scale = new Float32Array(traceInfo.npoints);
    scale.fill(pointSize);
    scatterPlot.setPointScaleFactors(scale);

    showLabels = showLabels && traceInfo.isCategorical;
    const labelsVisualizer = getVisualizer(scatterPlot, LABELS_VISUALIZER_ID);
    if (labelsVisualizer) {
        labelsVisualizer.labelsActive = showLabels;
        if (showLabels) {
            const labelsPositions = getCategoryLabelsPositions(traceInfo, categoricalNames);
            labelsVisualizer.fillStyle = darkMode ? 'white' : 'black';
            labelsVisualizer.shadowColor = !darkMode ? 'white' : 'black';
            labelsVisualizer.setLabels(labelsPositions.labels, labelsPositions.positions);
        }
    }

    scatterPlot.render();
}
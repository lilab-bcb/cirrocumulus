import {color} from 'd3-color';
import {makeStyles, ScatterPlot, ScatterPlotVisualizerSprites, ScatterPlotVisualizerSvgLabels} from 'scatter-gl';
import {Color, OrthographicCamera, Vector3} from 'three';
import {getVisualizer} from './ScatterChartThree';
import {getRgbScale, indexSort, randomSeq, rankIndexArray} from './util';

export const POINT_VISUALIZER_ID = 'SPRITES';

export const LABELS_VISUALIZER_ID = 'SVG_LABELS';

function scaleLinear(value, domain, range) {
    const domainDifference = domain[1] - domain[0];
    const rangeDifference = range[1] - range[0];

    const percentDomain = (value - domain[0]) / domainDifference;
    return percentDomain * rangeDifference + range[0];
}


export function getScaleFactor(size) {
    const ORTHO_CAMERA_FRUSTUM_HALF_EXTENT = 1.2;
    const aspectRatio = size.width / size.height;
    let left = -ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
    let right = ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
    let bottom = -ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
    let top = ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
    // Scale up the larger of (w, h) to match the aspect ratio.
    if (aspectRatio > 1) {
        left *= aspectRatio;
        right *= aspectRatio;
    } else {
        top /= aspectRatio;
        bottom /= aspectRatio;
    }
    let camera = new OrthographicCamera(
        left,
        right,
        top,
        bottom,
        -1000,
        1000
    );
    camera.up = new Vector3(0, 0, 1);

    camera.updateProjectionMatrix();
    return camera.projectionMatrix.elements[0];
}

export function createScatterPlot(containerElement, premultipliedAlpha, labels, interactive = true) {
    const styles = makeStyles();
    styles.label3D.fontSize = 40;

    const scatterPlot = new ScatterPlot(containerElement, {
        camera: {},
        selectEnabled: false,
        styles: styles,
        interactive: interactive

    }, premultipliedAlpha); // toDataUrl images are flipped on Safari when premultipliedAlpha is false
    let visualizers = [new ScatterPlotVisualizerSprites(styles)];
    if (labels) {
        // visualizers.push(new ScatterPlotVisualizerCanvasLabels(containerElement, styles));
        visualizers.push(new ScatterPlotVisualizerSvgLabels(containerElement, styles));
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
    let ranks;
    if (!is3d) {
        ranks = !trace.isCategorical ? rankIndexArray(indexSort(trace.values, true)) : randomSeq(trace.values.length, 1);
        // ranks go from 1 to values.length. Higher rank means higher value.
        zExtent[0] = 0;
        zExtent[1] = 1;
    }
    // Determine max and min of each axis of our data.
    for (let i = 0; i < npoints; i++) {
        const x = trace.x[i];
        if (x < xExtent[0]) {
            xExtent[0] = x;
        }
        if (x > xExtent[1]) {
            xExtent[1] = x;
        }

        const y = trace.y[i];
        if (y < yExtent[0]) {
            yExtent[0] = y;
        }
        if (y > yExtent[1]) {
            yExtent[1] = y;
        }
        if (is3d) {
            const z = trace.z[i];
            if (z < zExtent[0]) {
                zExtent[0] = z;
            }
            if (z > zExtent[1]) {
                zExtent[1] = z;
            }
        }
    }

    const getRange = (extent) => Math.abs(extent[1] - extent[0]);
    const xRange = getRange(xExtent);
    const yRange = getRange(yExtent);
    const zRange = getRange(zExtent);
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
        positions[dst++] = scaleLinear(is3d ? trace.z[i] : ranks[i] / (ranks.length + 1), zExtent, zScale);
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
            p = {x: [], y: [], z: []};
            categoryToPosition[value] = p;
            ncategories++;
        }
        p.count++;
        if (isImage) {
            p.x.push(traceInfo.x[i]);
            p.y.push(traceInfo.y[i]);
        } else {
            p.x.push(traceInfo.positions[j]);
            p.y.push(traceInfo.positions[j + 1]);
            p.z.push(traceInfo.positions[j + 2]);
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
        p.x.sort((a, b) => a - b);
        p.y.sort((a, b) => a - b);
        p.z.sort((a, b) => a - b);
        const mid = p.x.length / 2;
        labelPositions[positionIndex] = mid % 1 ? p.x[mid - 0.5] : (p.x[mid - 1] + p.x[mid]) / 2;
        labelPositions[positionIndex + 1] = mid % 1 ? p.y[mid - 0.5] : (p.y[mid - 1] + p.y[mid]) / 2;
        labelPositions[positionIndex + 2] = mid % 1 ? p.z[mid - 0.5] : (p.z[mid - 1] + p.z[mid]) / 2;
        positionIndex += 3;
    }

    return {labels: labelStrings, positions: labelPositions};
}

export function updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize, categoricalNames = {}, chartOptions) {
    const colors = traceInfo.colors;
    let positions = traceInfo.positions;

    const npoints = traceInfo.npoints;
    const isSelectionEmpty = selection.size === 0;

    for (let i = 0, j = 3, k = 2; i < npoints; i++, j += 4, k += 3) {
        const isSelected = isSelectionEmpty || selection.has(i);
        colors[j] = isSelected ? markerOpacity : unselectedMarkerOpacity;
    }
    scatterPlot.scene.background = chartOptions.darkMode ? new Color("rgb(0, 0, 0)") : null;
    scatterPlot.setDimensions(traceInfo.dimensions);
    let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
    spriteVisualizer.styles.fog.enabled = chartOptions.showFog;
    const axes = scatterPlot.scene.getObjectByName('axes');
    if (axes) {
        axes.visible = chartOptions.showAxis;
    }
    scatterPlot.setPointColors(colors);
    scatterPlot.setPointPositions(positions);

    // const {scaleDefault, scaleSelected, scaleHover} = scatterPlot.styles.point;

    const scale = new Float32Array(traceInfo.npoints);
    scale.fill(pointSize);
    scatterPlot.setPointScaleFactors(scale);

    const showLabels = chartOptions.showLabels && traceInfo.isCategorical;
    const labelsVisualizer = getVisualizer(scatterPlot, LABELS_VISUALIZER_ID);
    if (labelsVisualizer) {
        labelsVisualizer.labelsActive = showLabels;
        if (showLabels) {
            const labelsPositions = getCategoryLabelsPositions(traceInfo, categoricalNames);
            labelsVisualizer.fillStyle = chartOptions.darkMode ? 'white' : 'black';
            labelsVisualizer.shadowColor = chartOptions.darkMode ? 'rgba(0,0,0,0.9)' : 'rgba(255,255,255,0.9)';
            labelsVisualizer.shadowStroke = chartOptions.labelStrokeWidth;
            labelsVisualizer.setLabels(labelsPositions.labels, labelsPositions.positions);
            labelsVisualizer.font = 'bold ' + chartOptions.labelFontSize + 'px Roboto Condensed';
        }
    }

    scatterPlot.render();
}

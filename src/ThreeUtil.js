import {makeStyles, ScatterPlot, ScatterPlotVisualizerSprites, ScatterPlotVisualizerSvgLabels} from './scatter-gl';
import {Color, OrthographicCamera, Vector3} from 'three';
import {getEmbeddingKey} from './actions';
import {getVisualizer, setAxesColors} from './ScatterChartThree';
import {indexSort, randomSeq, rankIndexArray} from './util';

export const POINT_VISUALIZER_ID = 'SPRITES';
const SCATTER_PLOT_CUBE_LENGTH = 2;
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


export function getPositions(trace) {
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
        base * (range / maxRange)
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

export function getCategoryLabelsPositions(embedding, obsKeys, cachedData) {
    const obsArrayOfArrays = [];
    obsKeys.forEach(key => {
        obsArrayOfArrays.push(cachedData[key]);
    });
    const embeddingKey = getEmbeddingKey(embedding);
    const coordinates = cachedData[embeddingKey];
    if (coordinates == null) {
        throw new Error('Coordinates not found for ' + embedding.name);
    }
    const x = coordinates[embedding.name + '_1'];
    const y = coordinates[embedding.name + '_2'];
    const z = coordinates[embedding.name + '_3'];
    const is3d = z != null;
    const valueToCoords = {};
    let ncategories = 0;
    const isSpatial = embedding.spatial != null;
    const npoints = x.length;
    const nobs = obsKeys.length;
    for (let i = 0; i < npoints; i++) {
        let values = [];
        for (let j = 0; j < nobs; j++) {
            values.push(obsArrayOfArrays[j][i]);
        }

        const key = values.join(',');
        let p = valueToCoords[key];
        if (p === undefined) {
            p = {x: [], y: [], z: [], array: values};
            valueToCoords[key] = p;
            ncategories++;
        }
        p.count++;
        p.x.push(x[i]);
        p.y.push(y[i]);
        if (is3d) {
            p.z.push(z[i]);
        }
    }

    let xScale, yScale, zScale, xExtent, yExtent, zExtent;
    if (!isSpatial) {
        const getRange = (extent) => Math.abs(extent[1] - extent[0]);
        xExtent = [Infinity, -Infinity];
        yExtent = [Infinity, -Infinity];
        zExtent = is3d ? [Infinity, -Infinity] : [0, 1];
        // Determine max and min of each axis of our data.
        for (let i = 0; i < npoints; i++) {
            let value = x[i];
            if (value < xExtent[0]) {
                xExtent[0] = value;
            }
            if (value > xExtent[1]) {
                xExtent[1] = value;
            }

            value = y[i];
            if (value < yExtent[0]) {
                yExtent[0] = value;
            }
            if (value > yExtent[1]) {
                yExtent[1] = value;
            }
            if (is3d) {
                value = z[i];
                if (value < zExtent[0]) {
                    zExtent[0] = value;
                }
                if (value > zExtent[1]) {
                    zExtent[1] = value;
                }
            }
        }
        const xRange = getRange(xExtent);
        const yRange = getRange(yExtent);
        const zRange = getRange(zExtent);
        const maxRange = Math.max(xRange, yRange, zRange);
        const halfCube = SCATTER_PLOT_CUBE_LENGTH / 2;
        const makeScaleRange = (range, base) => [
            -base * (range / maxRange),
            base * (range / maxRange)
        ];
        xScale = makeScaleRange(xRange, halfCube);
        yScale = makeScaleRange(yRange, halfCube);
        zScale = makeScaleRange(zRange, halfCube);
    }

    let labelValues = [];
    let labelPositions = new Float32Array(ncategories * 3);
    let positionIndex = 0;
    for (let key in valueToCoords) {
        let p = valueToCoords[key];
        labelValues.push(p.array);
        p.x.sort((a, b) => a - b);
        p.y.sort((a, b) => a - b);
        p.z.sort((a, b) => a - b);
        const mid = p.x.length / 2;

        let xmedian = mid % 1 ? p.x[mid - 0.5] : (p.x[mid - 1] + p.x[mid]) / 2;
        let ymedian = mid % 1 ? p.y[mid - 0.5] : (p.y[mid - 1] + p.y[mid]) / 2;
        let zmedian = mid % 1 ? p.z[mid - 0.5] : (p.z[mid - 1] + p.z[mid]) / 2;

        if (!isSpatial) {
            xmedian = scaleLinear(xmedian, xExtent, xScale);
            ymedian = scaleLinear(ymedian, yExtent, yScale);
            zmedian = scaleLinear(zmedian, zExtent, zScale);
        }
        labelPositions[positionIndex] = xmedian;
        labelPositions[positionIndex + 1] = ymedian;
        labelPositions[positionIndex + 2] = is3d ? zmedian : 0;
        positionIndex += 3;
    }
    return {labels: labelValues, positions: labelPositions};
}

export function updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize, unselectedPointSize, categoricalNames = {}, chartOptions, obsCatKeys, cachedData, cameraDef) {
    const is3d = traceInfo.z != null;
    const colors = traceInfo.colors;
    let positions = traceInfo.positions;

    const npoints = traceInfo.x.length;
    const isSelectionEmpty = selection == null;
    const updateZ = !isSelectionEmpty && !is3d;
    if (updateZ) {
        positions = positions.slice();
    }
    const scale = new Float32Array(traceInfo.x.length);
    scale.fill(pointSize);
    for (let i = 0, j = 3, k = 2; i < npoints; i++, j += 4, k += 3) {
        const isSelected = isSelectionEmpty || selection.has(i);
        colors[j] = isSelected ? markerOpacity : unselectedMarkerOpacity;
        if (!isSelected) {
            scale[i] = unselectedPointSize;
        }
        if (updateZ && !isSelected) {
            positions[k] = -1;
        }
    }
    scatterPlot.scene.background = chartOptions.darkMode ? new Color("rgb(0, 0, 0)") : null;
    scatterPlot.setDimensions(traceInfo.dimensions);
    if (cameraDef) {
        scatterPlot.updateFromCameraDef(cameraDef);
    }
    let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
    spriteVisualizer.styles.fog.enabled = chartOptions.showFog;
    const axes = scatterPlot.scene.getObjectByName('axes');
    if (axes) {
        setAxesColors(scatterPlot, chartOptions.darkMode);
        axes.visible = chartOptions.showAxis;
    }
    scatterPlot.setPointColors(colors);
    scatterPlot.setPointPositions(positions);
    scatterPlot.setPointScaleFactors(scale);

    const showLabels = obsCatKeys.length > 0;
    const labelsVisualizer = getVisualizer(scatterPlot, LABELS_VISUALIZER_ID);
    if (labelsVisualizer) {
        labelsVisualizer.labelsActive = showLabels;
        if (showLabels) {
            const labelKey = getEmbeddingKey(traceInfo.embedding) + '_' + obsCatKeys.join(',');
            let labelsPositions = cachedData[labelKey];
            if (labelsPositions == null) {
                labelsPositions = getCategoryLabelsPositions(traceInfo.embedding, obsCatKeys, cachedData);
                cachedData[labelKey] = labelsPositions;
            }
            labelsVisualizer.fillStyle = chartOptions.darkMode ? 'white' : 'black';
            labelsVisualizer.shadowColor = chartOptions.darkMode ? 'rgba(0,0,0,0.9)' : 'rgba(255,255,255,0.9)';
            labelsVisualizer.shadowStroke = chartOptions.labelStrokeWidth;
            labelsVisualizer.font = 'bold ' + chartOptions.labelFontSize + 'px Roboto Condensed';
            const labels = getLabels(obsCatKeys, labelsPositions.labels, categoricalNames);

            labelsVisualizer.setLabels(labels, labelsPositions.positions);
        }
    }

    scatterPlot.render();
}

export function getLabels(obsCat, labels, categoricalNames) {
    let labelStrings = [];
    let renamedCategories = [];
    obsCat.forEach(key => renamedCategories.push(categoricalNames[key] || {}));
    for (let i = 0; i < labels.length; i++) {
        let array = labels[i];
        let value = [];
        for (let j = 0; j < array.length; j++) {
            const renamedValue = renamedCategories[j][array[j]];
            value.push(renamedValue != null && renamedValue.newValue != null ? renamedValue.newValue : array[j]);
        }
        labelStrings.push(value.join(','));
    }
    return labelStrings;
}

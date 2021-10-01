import withStyles from '@mui/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import React, {useEffect, useRef, useState} from 'react';
import {Color, Vector3, Vector4} from 'three';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {saveImage} from './ChartUtil';
import {numberFormat2f} from './formatters';
import {
    createScatterPlot,
    getCategoryLabelsPositions,
    getLabels,
    POINT_VISUALIZER_ID,
    updateScatterChart
} from './ThreeUtil';
import {indexSort, isPointInside} from './util';

function clamp(x, min_v, max_v) {
    return Math.min(Math.max(x, min_v), max_v);
}

function smoothstep(edge0, edge1, x) {
    const t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3 - 2 * t);
}

function mix(x, y, a) {
    return x * (1.0 - a) + y * a;
}

export function drawLabels(context, labels, positions, chartOptions, chartSize, camera) {
    const pos = new Vector3();
    context.fillStyle = chartOptions.darkMode ? 'white' : 'black';
    context.strokeStyle = chartOptions.darkMode ? 'rgba(0,0,0,0.9)' : 'rgba(255,255,255,0.9)';
    context.lineWidth = chartOptions.labelStrokeWidth;
    context.textAlign = 'center';
    context.textBaseline = "middle";
    const width = chartSize.width;
    const height = chartSize.height;
    const widthHalf = width / 2;
    const heightHalf = height / 2;
    for (let i = 0, k = 0; i < labels.length; i++, k += 3) {
        pos.x = positions[k];
        pos.y = positions[k + 1];
        pos.z = positions[k + 2];
        pos.project(camera);
        pos.x = (pos.x * widthHalf) + widthHalf;
        pos.y = -(pos.y * heightHalf) + heightHalf;

        context.strokeText(labels[i], pos.x, pos.y);
        context.fillText(labels[i], pos.x, pos.y);
    }
}

export function getVisualizer(scatterPlot, id) {
    for (let i = 0; i < scatterPlot.visualizers.length; i++) {
        if (scatterPlot.visualizers[i].id === id) {
            return scatterPlot.visualizers[i];
        }
    }
}

const styles = theme => ({

    root: {
        '& > *': {
            margin: theme.spacing(.4)
        },
        '& > .MuiIconButton-root': {
            padding: 0
        },
        '& > .cirro-active': {
            fill: 'rgb(220, 0, 78)',
            color: 'rgb(220, 0, 78)'
        },
        '& > .cirro-inactive': {
            fill: 'rgba(0, 0, 0, 0.26)',
            color: 'rgba(0, 0, 0, 0.26)'
        },
        position: 'absolute',
        top: 0,
        left: 0,
        display: 'inline-block',
        verticalAlign: 'top',
        whiteSpace: 'nowrap',
        overflow: 'hidden'
    }
});

export function setAxesColors(scatterPlot, darkMode) {
    const axes = scatterPlot.scene.getObjectByName('axes');
    if (axes) {
        axes.setColors(darkMode ? new Color("rgb(255, 255, 255)") : new Color("rgb(0, 0, 0)"));
    }
}

function ScatterChartThree(props) {
    const containerElementRef = useRef();
    const scatterPlotRef = useRef();
    const lastHoverIndexRef = useRef();
    const [forceUpdate, setForceUpdate] = useState(false);
    const previousChartSizeRef = useRef();

    function getSelectedIndex(point) {
        const trace = props.trace;
        const positions = trace.positions;
        const camera = scatterPlotRef.current.camera;
        const widthHalf = props.chartSize.width / 2;
        const heightHalf = props.chartSize.height / 2;
        const pos = new Vector3();
        let selectedIndex = -1;
        const tolerance = 2;
        if (lastHoverIndexRef.current !== -1) {
            pos.x = positions[lastHoverIndexRef.current * 3];
            pos.y = positions[lastHoverIndexRef.current * 3 + 1];
            pos.z = positions[lastHoverIndexRef.current * 3 + 2];
            pos.project(camera);
            pos.x = (pos.x * widthHalf) + widthHalf;
            pos.y = -(pos.y * heightHalf) + heightHalf;
            if (Math.abs(pos.x - point.x) <= tolerance && Math.abs(pos.y - point.y) <= tolerance) {
                selectedIndex = lastHoverIndexRef.current;
            }
        }

        if (selectedIndex === -1) {
            // TODO get all hover points
            for (let i = 0, k = 0, npoints = trace.x.length; i < npoints; i++, k += 3) {
                pos.x = positions[k];
                pos.y = positions[k + 1];
                pos.z = positions[k + 2];
                pos.project(camera);
                pos.x = (pos.x * widthHalf) + widthHalf;
                pos.y = -(pos.y * heightHalf) + heightHalf;
                if (Math.abs(pos.x - point.x) <= tolerance && Math.abs(pos.y - point.y) <= tolerance) {
                    selectedIndex = i;
                    break;
                }
            }
        }
        lastHoverIndexRef.current = selectedIndex;
        return selectedIndex;
    }

    useEffect(() => {
        init();
        if (previousChartSizeRef.current !== props.chartSize) {
            scatterPlotRef.current.resize();
        }
        draw();
        setAxesColors(scatterPlotRef.current, props.chartOptions.darkMode);
        props.chartOptions.scatterPlot = scatterPlotRef.current;
        if (props.chartOptions.camera) {
            scatterPlotRef.current.updateFromCameraDef(props.chartOptions.camera);
            props.chartOptions.camera = null;
        }
        previousChartSizeRef.current = props.chartSize;
        scatterPlotRef.current.clickCallback = (point, append) => {
            if (!props.trace.continuous) {
                const selectedIndex = getSelectedIndex(point);
                if (selectedIndex !== -1) {
                    props.handleClick({
                        name: props.trace.name,
                        value: props.trace.values[selectedIndex],
                        shiftKey: false,
                        metaKey: append
                    });
                }
                lastHoverIndexRef.current = selectedIndex;
            }

        };
        scatterPlotRef.current.hoverCallback = (point) => {
            if (point == null) {
                props.setTooltip('');
            } else {
                const selectedIndex = getSelectedIndex(point);
                if (selectedIndex !== -1) {
                    const {trace} = props;
                    let value = trace.values[selectedIndex];
                    let categoryObject = props.categoricalNames[trace.name] || {};
                    let renamedValue = categoryObject[value];
                    if (renamedValue != null && renamedValue.newValue != null) {
                        value = renamedValue.newValue;
                    }

                    if (typeof value === 'number') {
                        value = numberFormat2f(value);
                        if (value.endsWith('.00')) {
                            value = value.substring(0, value.lastIndexOf('.'));
                        }
                    }
                    props.setTooltip('' + value);
                } else {
                    props.setTooltip('');
                }
            }
        };
        scatterPlotRef.current.lassoCallback = (points, appendToSelection) => {
            const trace = props.trace;
            const positions = trace.positions;
            const camera = scatterPlotRef.current.camera;
            const widthHalf = props.chartSize.width / 2;
            const heightHalf = props.chartSize.height / 2;
            const pos = new Vector3();
            const selectedIndices = new Set();

            for (let i = 0, k = 0, npoints = trace.x.length; i < npoints; i++, k += 3) {
                pos.x = positions[k];
                pos.y = positions[k + 1];
                pos.z = positions[k + 2];
                pos.project(camera);
                pos.x = (pos.x * widthHalf) + widthHalf;
                pos.y = -(pos.y * heightHalf) + heightHalf;
                if (isPointInside(pos, points)) {
                    selectedIndices.add(i);
                }
            }
            if (selectedIndices.size === 0) {
                props.onSelected({name: getEmbeddingKey(trace.embedding)});
            } else {
                props.onSelected({
                    name: getEmbeddingKey(trace.embedding),
                    clear: !appendToSelection,
                    value: {basis: trace.embedding, indices: selectedIndices}
                });
            }
        };
        scatterPlotRef.current.boxCallback = (rect, appendToSelection) => {
            if (scatterPlotRef.current.interactionMode === 'PAN') {
                return;
            }
            const trace = props.trace;
            const positions = trace.positions;
            const camera = scatterPlotRef.current.camera;
            const widthHalf = props.chartSize.width / 2;
            const heightHalf = props.chartSize.height / 2;
            const pos = new Vector3();
            const selectedIndices = new Set();

            for (let i = 0, k = 0, npoints = trace.x.length; i < npoints; i++, k += 3) {
                pos.x = positions[k];
                pos.y = positions[k + 1];
                pos.z = positions[k + 2];
                pos.project(camera);
                pos.x = (pos.x * widthHalf) + widthHalf;
                pos.y = -(pos.y * heightHalf) + heightHalf;
                if (pos.x >= rect.x && pos.x <= (rect.x + rect.width) && pos.y >= rect.y && pos.y <= (rect.y + rect.height)) {
                    selectedIndices.add(i);
                }
            }

            if (selectedIndices.size === 0) {
                props.onSelected({name: getEmbeddingKey(trace.embedding)});
            } else {
                props.onSelected({
                    name: getEmbeddingKey(trace.embedding),
                    clear: !appendToSelection,
                    value: {basis: trace.embedding, indices: selectedIndices}
                });
            }
        };
        scatterPlotRef.current.cameraCallback = (eventName) => {
            if (scatterPlotRef.current.interactionMode === 'PAN' && props.trace.dimensions === 3) {
                // repaint gallery charts with same embedding
                if (eventName === 'end') {
                    cameraCallback(eventName);
                }
            }
        };
    });

    useEffect(() => {
        return () => {
            scatterPlotRef.current.dispose();
        };
    }, []);

    function calculatePointSize(trace) {
        const n = trace.x.length;
        const SCALE = 200;
        const LOG_BASE = 8;
        const DIVISOR = 1.5;
        // Scale point size inverse-logarithmically to the number of points.
        const pointSize = SCALE / Math.log(n) / Math.log(LOG_BASE);
        return trace.dimensions === 3 ? pointSize : pointSize / DIVISOR;
    }

    function drawContext(context, chartSize, format) {
        const {
            obsCat,
            cachedData,
            trace,
            markerOpacity,
            unselectedMarkerOpacity,
            selection,
            categoricalNames,
            chartOptions
        } = props;
        const scatterPlot = scatterPlotRef.current;
        const pointSize = calculatePointSize(trace);
        const scaleFactor = props.pointSize;
        const PI2 = 2 * Math.PI;
        const colors = trace.colors;
        const positions = trace.positions;
        const camera = scatterPlot.camera;
        const width = chartSize.width;
        const height = chartSize.height;
        if (chartOptions.darkMode) {
            context.fillStyle = 'black';
            context.fillRect(0, 0, width, height);
        }
        const widthHalf = width / 2;
        const heightHalf = height / 2;
        const colorScaleConverter = scaleLinear().domain([0, 1]).range([0, 255]);
        const npoints = trace.x.length;
        const is3d = trace.dimensions === 3;
        let outputPointSize;
        let fog = scatterPlot.scene.fog;
        let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
        // const zoomFactor = getScaleFactor(props.chartSize);
        // const zoomFactorSpecified = false;

        if (!is3d) {
            const PI = 3.1415926535897932384626433832795;
            const minScale = 0.1;  // minimum scaling factor
            const outSpeed = 2.0;  // shrink speed when zooming out
            const outNorm = (1. - minScale) / Math.atan(outSpeed);
            const maxScale = 15.0;  // maximum scaling factor
            const inSpeed = 0.02;  // enlarge speed when zooming in
            const zoomOffset = 0.3;  // offset zoom pivot
            let m = camera.projectionMatrix.elements[0];
            // if (zoomFactorSpecified) {
            //     m = zoomFactor;
            // }
            let zoom = m + zoomOffset;  // zoom pivot
            let scale = zoom < 1. ? 1. + outNorm * Math.atan(outSpeed * (zoom - 1.)) :
                1. + 2. / PI * (maxScale - 1.) * Math.atan(inSpeed * (zoom - 1.));
            outputPointSize = pointSize * scale;
        }
        let gl_PointSize = (outputPointSize * scaleFactor);
        gl_PointSize /= 2;
        const pos = new Vector3();
        let cameraSpacePos = new Vector4();
        let object = spriteVisualizer.points;
        let modelViewMatrix = object.modelViewMatrix.clone();
        modelViewMatrix.multiplyMatrices(camera.matrixWorldInverse, object.matrixWorld);
        const showFog = chartOptions.showFog;
        const isSelectionEmpty = selection.size === 0;
        let pointOrder;
        if (is3d) {
            pointOrder = new Uint32Array(npoints);
            for (let i = 0; i < npoints; i++) {
                pointOrder[i] = i;
            }
        } else {
            const z = new Float32Array(npoints);
            for (let i = 0; i < npoints; i++) {
                z[i] = positions[i * 3 + 2];
            }
            pointOrder = indexSort(z, true);
        }
        for (let j = 0; j < npoints; j++) {
            const index = pointOrder[j];
            const positionIndex = index * 3;
            const colorIndex = index * 4;
            const isSelected = isSelectionEmpty || selection.has(index);
            pos.x = positions[positionIndex];
            pos.y = positions[positionIndex + 1];
            pos.z = positions[positionIndex + 2];
            pos.project(camera);

            let r = colors[colorIndex];
            let g = colors[colorIndex + 1];
            let b = colors[colorIndex + 2];
            let a = isSelected ? markerOpacity : unselectedMarkerOpacity;
            if (is3d) {
                cameraSpacePos.x = positions[positionIndex];
                cameraSpacePos.y = positions[positionIndex + 1];
                cameraSpacePos.z = positions[positionIndex + 2];
                cameraSpacePos.w = 1;
                cameraSpacePos.applyMatrix4(modelViewMatrix);
                outputPointSize = -pointSize / cameraSpacePos.z;
                gl_PointSize = (outputPointSize * scaleFactor) / 4;
                if (showFog) {
                    const fogDepth = pointSize / outputPointSize * 1.2;
                    const fogFactor = smoothstep(fog.near, fog.far, fogDepth);
                    r = mix(r, fog.color.r, fogFactor);
                    g = mix(g, fog.color.g, fogFactor);
                    b = mix(b, fog.color.b, fogFactor);
                }
            }
            pos.x = (pos.x * widthHalf) + widthHalf;
            pos.y = -(pos.y * heightHalf) + heightHalf;

            r = Math.round(colorScaleConverter(r));
            g = Math.round(colorScaleConverter(g));
            b = Math.round(colorScaleConverter(b));

            context.fillStyle = 'rgba(' + r + ',' + g + ',' + b + ',' + a + ')';
            context.beginPath();
            context.arc(pos.x, pos.y, gl_PointSize, 0, PI2);
            context.closePath();
            context.fill();
        }
        if (obsCat.length > 0) {
            const labelsPositions = getCategoryLabelsPositions(trace.embedding, obsCat, cachedData);
            let font = format === 'svg' ? 'serif' : 'Roboto Condensed';
            context.font = 'bold ' + chartOptions.labelFontSize + 'px ' + font;
            drawLabels(context, getLabels(obsCat, labelsPositions.labels, categoricalNames), labelsPositions.positions, chartOptions, chartSize, camera);
        }
    }

    function onSaveImage(format) {
        const {trace, chartSize} = props;
        saveImage(trace, chartSize, drawContext, format);
    }


    function resetCamera() {
        const scatterPlot = scatterPlotRef.current;
        scatterPlot.resetZoom();
        if (scatterPlot.interactionMode === 'PAN' && props.trace.dimensions === 3) {
            props.onCamera('change', scatterPlot.getCameraDef());
        }
    }


    function onShowAxis() {
        const scatterPlot = scatterPlotRef.current;
        const axes = scatterPlot.scene.getObjectByName('axes');
        props.chartOptions.showAxis = !props.chartOptions.showAxis;
        if (axes) {
            axes.visible = props.chartOptions.showAxis;
        }

        props.setChartOptions(props.chartOptions);
    }

    function onShowFog() {
        const scatterPlot = scatterPlotRef.current;
        props.chartOptions.showFog = !props.chartOptions.showFog;
        const spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
        spriteVisualizer.styles.fog.enabled = props.chartOptions.showFog;
        props.setChartOptions(props.chartOptions);
    }


    function onDragMode(mode) {
        const scatterPlot = scatterPlotRef.current;
        if (mode === 'pan') {
            scatterPlot.setInteractionMode('PAN');
        } else if (mode === 'select') {
            scatterPlot.setInteractionMode('SELECT');
            scatterPlot.rectangleSelector.setSelectionMode('BOX');
        } else if (mode === 'lasso') {
            scatterPlot.rectangleSelector.setSelectionMode('LASSO');
            scatterPlot.setInteractionMode('SELECT');
        }
        props.chartOptions.dragmode = mode;
        props.setChartOptions(props.chartOptions);
    }

    function onToggleAnimation() {
        const scatterPlot = scatterPlotRef.current;
        if (scatterPlot.orbitIsAnimating()) {
            scatterPlot.stopOrbitAnimation();
            props.chartOptions.animating = false;
        } else {
            scatterPlot.startOrbitAnimation();
            props.chartOptions.animating = true;
        }
        props.setChartOptions(props.chartOptions);

    }

    function cameraCallback(eventName) {
        props.onCamera(eventName, scatterPlotRef.current.getCameraDef());
    }

    function init() {
        if (scatterPlotRef.current == null) {
            scatterPlotRef.current = createScatterPlot(containerElementRef.current, window.ApplePaySession, true);
            if (props.chartOptions.dragmode === 'pan') {
                scatterPlotRef.current.setInteractionMode('PAN');
            } else if (props.chartOptions.dragmode === 'select') {
                scatterPlotRef.current.setInteractionMode('SELECT');
                scatterPlotRef.current.rectangleSelector.setSelectionMode('BOX');
            } else if (props.chartOptions.dragmode === 'lasso') {
                scatterPlotRef.current.setInteractionMode('SELECT');
                scatterPlotRef.current.rectangleSelector.setSelectionMode('LASSO');
            }
            const axes = scatterPlotRef.current.scene.getObjectByName('axes');
            if (axes) {
                axes.visible = props.chartOptions.showAxis;
                axes.setColors(props.chartOptions.darkMode ? new Color("rgb(255, 255, 255)") : new Color("rgb(0, 0, 0)"));
            }
            getVisualizer(scatterPlotRef.current, POINT_VISUALIZER_ID).styles.fog.enabled = props.chartOptions.showFog;
            const canvas = containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
            const webglcontextlost = (e) => {
                console.log('lost webgl context');
                e.preventDefault();
                scatterPlotRef.current = null;
            };
            const webglcontextrestored = (e) => {
                console.log('restored webgl context');
                e.preventDefault();
                setForceUpdate(!forceUpdate);
            };
            canvas.addEventListener('webglcontextlost', webglcontextlost);
            canvas.addEventListener('webglcontextrestored', webglcontextrestored);
            return true;
        }
        return false;

    }


    function draw() {
        const {
            obsCat,
            cachedData,
            trace,
            markerOpacity,
            unselectedMarkerOpacity,
            selection,
            pointSize,
            chartOptions,
            categoricalNames,
            unselectedPointSize
        } = props;
        updateScatterChart(scatterPlotRef.current, trace, selection, markerOpacity, unselectedMarkerOpacity, pointSize, unselectedPointSize,
            categoricalNames, chartOptions, obsCat, cachedData);
    }


    return <>
        <div className={props.classes.root}>
            <ChartToolbar
                dragmode={props.chartOptions.dragmode}
                // editSelection={props.chartOptions.editSelection}
                onGallery={props.onGallery}
                animating={props.chartOptions.animating}
                showFog={props.chartOptions.showFog}
                onShowFog={onShowFog}
                is3d={props.trace && props.trace.z != null}
                toggleAnimation={onToggleAnimation}
                onSaveImage={onSaveImage}
                onDragMode={onDragMode}
                // onCopyImage={onCopyImage}
                // onEditSelection={onEditSelection}
                onShowAxis={onShowAxis}
                onHome={resetCamera}
                showAxis={props.chartOptions.showAxis}/>
        </div>
        <div data-testid="scatter-chart-three" style={{
            display: 'inline-block',
            width: props.chartSize.width,
            height: props.chartSize.height
        }}
             ref={containerElementRef}>
        </div>
    </>;

}

export default withStyles(styles)(ScatterChartThree);



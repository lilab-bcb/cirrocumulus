import withStyles from '@mui/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import React, {useEffect, useRef, useState} from 'react';
import {Color, Vector3, Vector4} from 'three';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {numberFormat2f} from './formatters';
import {
    createScatterPlot,
    getCategoryLabelsPositions,
    getLabels,
    POINT_VISUALIZER_ID,
    updateScatterChart
} from './ThreeUtil';
import {indexSort, isPointInside} from './util';
import {saveImage} from './ChartUtil';
import CirroTooltip, {SCATTER_TRANSITION} from './CirroTooltip';

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
    const previousCameraPosition = useRef({x: -1, y: -1, z: -1});
    const [tip, setTip] = useState({html: ''});

    const {
        cachedData,
        categoricalNames,
        chartOptions,
        chartSize,
        classes,
        handleClick,
        markerOpacity,
        obsCat,
        onCamera,
        onGallery,
        onSelected,
        pointSize,
        selection,
        setChartOptions,
        trace,
        unselectedMarkerOpacity,
        unselectedPointSize
    } = props;

    const darkMode = chartOptions.darkMode;
    const showAxis = chartOptions.showAxis;
    const showFog = chartOptions.showFog;

    useEffect(() => {
        if (scatterPlotRef.current == null) {
            const dragmode = chartOptions.dragmode;
            scatterPlotRef.current = createScatterPlot(containerElementRef.current, window.ApplePaySession, true);
            if (dragmode === 'pan') {
                scatterPlotRef.current.setInteractionMode('PAN');
            } else if (dragmode === 'select') {
                scatterPlotRef.current.setInteractionMode('SELECT');
                scatterPlotRef.current.rectangleSelector.setSelectionMode('BOX');
            } else if (dragmode === 'lasso') {
                scatterPlotRef.current.setInteractionMode('SELECT');
                scatterPlotRef.current.rectangleSelector.setSelectionMode('LASSO');
            }

            const canvas = containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
            const webglcontextlost = (e) => {
                console.log('lost webgl context');
                e.preventDefault();
            };
            const webglcontextrestored = (e) => {
                console.log('restored webgl context');
                e.preventDefault();
                setForceUpdate(c => !c);
            };
            canvas.addEventListener('webglcontextlost', webglcontextlost);
            canvas.addEventListener('webglcontextrestored', webglcontextrestored);
        }
        if (chartOptions.camera) {
            scatterPlotRef.current.updateFromCameraDef(chartOptions.camera);
            chartOptions.camera = null;
        }
        chartOptions.scatterPlot = scatterPlotRef.current;
    }, [scatterPlotRef, containerElementRef, chartOptions]);


    useEffect(() => {

        function handleCamera(eventName) {
            onCamera(eventName, scatterPlotRef.current.getCameraDef());
        }

        function getSelectedIndex(point) {
            const positions = trace.positions;
            const camera = scatterPlotRef.current.camera;
            const widthHalf = chartSize.width / 2;
            const heightHalf = chartSize.height / 2;
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


        let singleClickTimer;
        let clickCount = 0;
        scatterPlotRef.current.clickCallback = (point, append) => {
            if (!trace.continuous) {
                clickCount++;
                if (clickCount === 1) {
                    singleClickTimer = setTimeout(() => {
                        clickCount = 0; // reset after timeout
                    }, 400);
                } else if (clickCount === 2) {
                    clearTimeout(singleClickTimer);
                    clickCount = 0;
                    const selectedIndex = getSelectedIndex(point);
                    if (selectedIndex !== -1) {
                        handleClick({
                            name: trace.name,
                            value: trace.values[selectedIndex],
                            shiftKey: false,
                            metaKey: append
                        });
                    }
                    lastHoverIndexRef.current = selectedIndex;
                }
            }
        };
        scatterPlotRef.current.hoverCallback = (point, event) => {
            if (point == null) {
                setTip({html: ''});
            } else {
                const selectedIndex = getSelectedIndex(point);
                if (selectedIndex !== -1) {
                    let value = trace.values[selectedIndex];
                    let categoryObject = categoricalNames[trace.name] || {};
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
                    setTip({html: '' + value, clientX: event.clientX, clientY: event.clientY});
                    // updateTooltipText(tip, '' + value, event);
                } else {
                    setTip({html: ''});
                    // updateTooltipText(tip, '', event);
                }
            }
        };
        scatterPlotRef.current.lassoCallback = (points, appendToSelection) => {
            const positions = trace.positions;
            const camera = scatterPlotRef.current.camera;
            const widthHalf = chartSize.width / 2;
            const heightHalf = chartSize.height / 2;
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
                onSelected({name: getEmbeddingKey(trace.embedding)});
            } else {
                onSelected({
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
            const positions = trace.positions;
            const camera = scatterPlotRef.current.camera;
            const widthHalf = chartSize.width / 2;
            const heightHalf = chartSize.height / 2;
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
                onSelected({name: getEmbeddingKey(trace.embedding)});
            } else {
                onSelected({
                    name: getEmbeddingKey(trace.embedding),
                    clear: !appendToSelection,
                    value: {basis: trace.embedding, indices: selectedIndices}
                });
            }
        };
        scatterPlotRef.current.cameraCallback = (eventName, position, target) => {
            if (scatterPlotRef.current.interactionMode === 'PAN' && trace.dimensions === 3) {
                // repaint gallery charts with same embedding
                if (eventName === 'end' && (previousCameraPosition.current.x != position.x || previousCameraPosition.current.y != position.y || previousCameraPosition.current.z != position.z)) {
                    previousCameraPosition.current = {x: position.x, y: position.y, z: position.z};
                    handleCamera(eventName);
                }
            }
        };
    }, [scatterPlotRef, categoricalNames, chartSize, trace]); // onSelected, handleClick, onCamera

    useEffect(() => {
        setAxesColors(scatterPlotRef.current, darkMode);
        const axes = scatterPlotRef.current.scene.getObjectByName('axes');
        if (axes) {
            axes.visible = showAxis;
        }
        getVisualizer(scatterPlotRef.current, POINT_VISUALIZER_ID).styles.fog.enabled = showFog;
    }, [scatterPlotRef, darkMode, showAxis, showFog]);


    useEffect(() => {
        if (previousChartSizeRef.current !== chartSize) {
            scatterPlotRef.current.resize();
        }
        previousChartSizeRef.current = chartSize;
        updateScatterChart(scatterPlotRef.current, trace, selection, markerOpacity, unselectedMarkerOpacity, pointSize, unselectedPointSize,
            categoricalNames, chartOptions, obsCat, cachedData);
    }, [scatterPlotRef, trace, selection, markerOpacity, unselectedMarkerOpacity, pointSize, unselectedPointSize,
        categoricalNames, chartOptions, obsCat, cachedData, chartSize, forceUpdate]);

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
        const PI2 = 2 * Math.PI;
        const width = chartSize.width;
        const height = chartSize.height;
        const scatterPlot = scatterPlotRef.current;
        const colors = trace.colors;
        const positions = trace.positions;


        if (chartOptions.darkMode) {
            context.fillStyle = 'black';
            context.fillRect(0, 0, width, height);
        }
        const widthHalf = width / 2;
        const heightHalf = height / 2;
        const colorScaleConverter = scaleLinear().domain([0, 1]).range([0, 255]);
        const npoints = trace.x.length;
        const is3d = trace.dimensions === 3;
        const fog = scatterPlot.scene.fog;
        const camera = scatterPlot.camera;
        const spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
        const scaleFactor = is3d ? calculatePointSize(trace) : 1; // TODO
        const pos = new Vector3();
        const cameraSpacePos = new Vector4();
        const modelViewMatrix = spriteVisualizer.points.modelViewMatrix.clone();
        modelViewMatrix.multiplyMatrices(camera.matrixWorldInverse, spriteVisualizer.points.matrixWorld);
        const showFog = chartOptions.showFog && is3d;
        const isSelectionEmpty = selection == null;
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
            let outputPointSize;
            if (is3d) {
                cameraSpacePos.x = positions[positionIndex];
                cameraSpacePos.y = positions[positionIndex + 1];
                cameraSpacePos.z = positions[positionIndex + 2];
                cameraSpacePos.w = 1;
                cameraSpacePos.applyMatrix4(modelViewMatrix);
                outputPointSize = -pointSize / cameraSpacePos.z;
            } else {  // Create size attenuation (if we're in 2D mode)
                const PI = 3.1415926535897932384626433832795;
                const minScale = 0.1;  // minimum scaling factor
                const outSpeed = 2.0;  // shrink speed when zooming out
                const outNorm = (1. - minScale) / Math.atan(outSpeed);
                const maxScale = 15.0;  // maximum scaling factor
                const inSpeed = 0.02;  // enlarge speed when zooming in
                const zoomOffset = 0.3;  // offset zoom pivot

                cameraSpacePos.x = positions[positionIndex];
                cameraSpacePos.y = positions[positionIndex + 1];
                cameraSpacePos.z = positions[positionIndex + 2];
                cameraSpacePos.w = 1;
                cameraSpacePos.applyMatrix4(modelViewMatrix);
                let m = -cameraSpacePos.z;
                // if (zoomFactorSpecified) {
                //     m = zoomFactor;
                // }
                const zoom = m + zoomOffset;  // zoom pivot
                const scale = zoom < 1. ? 1. + outNorm * Math.atan(outSpeed * (zoom - 1.)) :
                    1. + 2. / PI * (maxScale - 1.) * Math.atan(inSpeed * (zoom - 1.));
                outputPointSize = pointSize * scale;
            }
            if (showFog) {
                const fogDepth = pointSize / outputPointSize * 1.2;
                const fogFactor = smoothstep(fog.near, fog.far, fogDepth);
                r = mix(r, fog.color.r, fogFactor);
                g = mix(g, fog.color.g, fogFactor);
                b = mix(b, fog.color.b, fogFactor);
            }
            const gl_PointSize = (outputPointSize * scaleFactor) / 2;
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
            const font = format === 'svg' ? 'serif' : 'Roboto Condensed';
            context.font = 'bold ' + chartOptions.labelFontSize + 'px ' + font;
            drawLabels(context, getLabels(obsCat, labelsPositions.labels, categoricalNames), labelsPositions.positions, chartOptions, chartSize, camera);
        }
    }

    function onSaveImage(format) {
        saveImage(trace, chartSize, drawContext, format);
    }


    function resetCamera() {
        const scatterPlot = scatterPlotRef.current;
        scatterPlot.resetZoom();
        if (scatterPlot.interactionMode === 'PAN' && trace.dimensions === 3) {
            onCamera('change', scatterPlot.getCameraDef());
        }
    }


    function onShowAxis() {
        const scatterPlot = scatterPlotRef.current;
        const axes = scatterPlot.scene.getObjectByName('axes');
        chartOptions.showAxis = !chartOptions.showAxis;
        if (axes) {
            axes.visible = chartOptions.showAxis;
        }

        setChartOptions(chartOptions);
    }

    function onShowFog() {
        const scatterPlot = scatterPlotRef.current;
        chartOptions.showFog = !chartOptions.showFog;
        const spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
        spriteVisualizer.styles.fog.enabled = chartOptions.showFog;
        setChartOptions(chartOptions);
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
        chartOptions.dragmode = mode;
        setChartOptions(chartOptions);
    }

    function onToggleAnimation() {
        const scatterPlot = scatterPlotRef.current;
        if (scatterPlot.orbitIsAnimating()) {
            scatterPlot.stopOrbitAnimation();
            chartOptions.animating = false;
        } else {
            scatterPlot.startOrbitAnimation();
            chartOptions.animating = true;
        }
        setChartOptions(chartOptions);

    }


    return <>
        <div className={classes.root}>
            <ChartToolbar
                dragmode={chartOptions.dragmode}
                // editSelection={chartOptions.editSelection}
                onGallery={onGallery}
                animating={chartOptions.animating}
                showFog={chartOptions.showFog}
                onShowFog={onShowFog}
                is3d={trace && trace.z != null}
                toggleAnimation={onToggleAnimation}
                onSaveImage={onSaveImage}
                onDragMode={onDragMode}
                // onCopyImage={onCopyImage}
                // onEditSelection={onEditSelection}
                onShowAxis={onShowAxis}
                onHome={resetCamera}
                showAxis={chartOptions.showAxis}/>
        </div>
        <div data-testid="scatter-chart-three" style={{
            display: 'inline-block',
            width: chartSize.width,
            height: chartSize.height
        }}
             ref={containerElementRef}>
            <CirroTooltip html={tip.html} clientX={tip.clientX} clientY={tip.clientY} transition={SCATTER_TRANSITION}/>
        </div>
    </>;

}

export default withStyles(styles)(ScatterChartThree);



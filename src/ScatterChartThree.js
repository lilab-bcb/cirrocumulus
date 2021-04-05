import withStyles from '@material-ui/core/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import {bind} from 'lodash';
import React from 'react';
import {Vector3, Vector4} from 'three';
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
};

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
            margin: theme.spacing(.4),
        },
        '& > .MuiIconButton-root': {
            padding: 0,
        },
        '& > .cirro-active': {
            fill: 'rgb(220, 0, 78)',
            color: 'rgb(220, 0, 78)',
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
        overflow: 'hidden',
    }
});

class ScatterChartThree extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();
        this.scatterPlot = null;
        this.lastHoverIndex = -1;
        this.state = {forceUpdate: false};
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.chartSize !== this.props.chartSize) {
            this.scatterPlot.resize();
        }
        this.init();
        this.draw();
    }

    componentDidMount() {
        this.init();
        this.draw(); // TODO fix 2x draw
        this.draw();
        if (this.props.chartOptions.camera) {
            this.scatterPlot.updateFromCameraDef(this.props.chartOptions.camera);
            this.props.chartOptions.camera = null;
            this.draw();
        }
    }


    componentWillUnmount() {
        // if (this.scatterGL != null) {
        //     this.scatterGL.dispose();
        // }
    }


    calculatePointSize(trace) {
        const n = trace.npoints;
        const SCALE = 200;
        const LOG_BASE = 8;
        const DIVISOR = 1.5;
        // Scale point size inverse-logarithmically to the number of points.
        const pointSize = SCALE / Math.log(n) / Math.log(LOG_BASE);
        return trace.dimensions === 3 ? pointSize : pointSize / DIVISOR;
    }

    drawContext(context, chartSize, format) {
        const {
            obsCat,
            cachedData,
            trace,
            markerOpacity,
            unselectedMarkerOpacity,
            selection,
            categoricalNames,
            chartOptions
        } = this.props;
        const pointSize = this.calculatePointSize(trace);
        const scaleFactor = this.props.pointSize;
        const PI2 = 2 * Math.PI;
        const colors = trace.colors;
        const positions = trace.positions;
        const camera = this.scatterPlot.camera;
        const width = chartSize.width;
        const height = chartSize.height;
        if (chartOptions.darkMode) {
            context.fillStyle = 'black';
            context.fillRect(0, 0, width, height);
        }
        const widthHalf = width / 2;
        const heightHalf = height / 2;
        const colorScale = scaleLinear().domain([0, 1]).range([0, 255]);
        const npoints = trace.npoints;
        const is3d = trace.dimensions === 3;
        let outputPointSize;
        let fog = this.scatterPlot.scene.fog;
        let spriteVisualizer = getVisualizer(this.scatterPlot, POINT_VISUALIZER_ID);

        if (!is3d) {
            const PI = 3.1415926535897932384626433832795;
            const minScale = 0.1;  // minimum scaling factor
            const outSpeed = 2.0;  // shrink speed when zooming out
            const outNorm = (1. - minScale) / Math.atan(outSpeed);
            const maxScale = 15.0;  // maximum scaling factor
            const inSpeed = 0.02;  // enlarge speed when zooming in
            const zoomOffset = 0.3;  // offset zoom pivot
            let zoom = camera.projectionMatrix.elements[0] + zoomOffset;  // zoom pivot
            let scale = zoom < 1. ? 1. + outNorm * Math.atan(outSpeed * (zoom - 1.)) :
                1. + 2. / PI * (maxScale - 1.) * Math.atan(inSpeed * (zoom - 1.));
            outputPointSize = pointSize * scale;
        }
        let gl_PointSize = (outputPointSize * scaleFactor) / 4;
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

            r = Math.round(colorScale(r));
            g = Math.round(colorScale(g));
            b = Math.round(colorScale(b));

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

    // onCopyImage = (event) => {
    //     const size = this.chartSize;
    //     const canvas = document.createElement('canvas');
    //     canvas.width = size.width; // * window.devicePixelRatio;
    //     canvas.height = size.height; // * window.devicePixelRatio;
    //     const context = canvas.getContext('2d');
    //     this.drawContext(context, this.chartSize);
    //     const url = canvas.toDataURL();
    //     setClipboardData([{
    //         format: 'image/png',
    //         data: '<img src="' + url + '">'
    //     }]);
    // };


    onSaveImage = (format) => {
        const {trace, chartSize} = this.props;
        saveImage(trace, chartSize, bind(this.drawContext, this), format);
    };


    // onEditSelection = () => {
    //     this.props.chartOptions.editSelection = !this.props.chartOptions.editSelection;
    //     this.props.setChartOptions(this.props.chartOptions);
    // };


    resetCamera = () => {
        this.scatterPlot.resetZoom();
    };


    onShowAxis = () => {
        const axes = this.scatterPlot.scene.getObjectByName('axes');
        this.props.chartOptions.showAxis = !this.props.chartOptions.showAxis;
        if (axes) {
            axes.visible = this.props.chartOptions.showAxis;
        }

        this.props.setChartOptions(this.props.chartOptions);
    };

    onShowFog = () => {
        let spriteVisualizer = getVisualizer(this.scatterPlot, POINT_VISUALIZER_ID);
        this.props.chartOptions.showFog = !this.props.chartOptions.showFog;
        spriteVisualizer.styles.fog.enabled = this.props.chartOptions.showFog;
        this.props.setChartOptions(this.props.chartOptions);
    };


    onDragMode = (mode) => {

        if (mode === 'pan') {
            this.scatterPlot.setInteractionMode('PAN');
        } else if (mode === 'select') {
            this.scatterPlot.setInteractionMode('SELECT');
            this.scatterPlot.rectangleSelector.setSelectionMode('BOX');
        } else if (mode === 'lasso') {
            this.scatterPlot.rectangleSelector.setSelectionMode('LASSO');
            this.scatterPlot.setInteractionMode('SELECT');
        }
        this.props.chartOptions.dragmode = mode;
        this.props.setChartOptions(this.props.chartOptions);
    };

    onToggleAnimation = () => {
        if (this.scatterPlot.orbitIsAnimating()) {
            this.scatterPlot.stopOrbitAnimation();
            this.props.chartOptions.animating = false;
        } else {
            this.scatterPlot.startOrbitAnimation();
            this.props.chartOptions.animating = true;
        }
        this.props.setChartOptions(this.props.chartOptions);
    };


    init() {
        if (this.scatterPlot == null) {
            const containerElement = this.containerElementRef.current;
            this.scatterPlot = createScatterPlot(containerElement, window.ApplePaySession, true);
            if (this.props.chartOptions.dragmode === 'pan') {
                this.scatterPlot.setInteractionMode('PAN');
            } else if (this.props.chartOptions.dragmode === 'select') {
                this.scatterPlot.setInteractionMode('SELECT');
                this.scatterPlot.rectangleSelector.setSelectionMode('BOX');
            } else if (this.props.chartOptions.dragmode === 'lasso') {
                this.scatterPlot.setInteractionMode('SELECT');
                this.scatterPlot.rectangleSelector.setSelectionMode('LASSO');

            }
            const axes = this.scatterPlot.scene.getObjectByName('axes');
            if (axes) {
                axes.visible = this.props.chartOptions.showAxis;
            }
            let spriteVisualizer = getVisualizer(this.scatterPlot, POINT_VISUALIZER_ID);
            spriteVisualizer.styles.fog.enabled = this.props.chartOptions.showFog;
            this.scatterPlot.hoverCallback = (point) => {
                if (point == null) {
                    this.props.setTooltip('');
                } else {
                    const trace = this.props.trace;
                    const positions = trace.positions;
                    const camera = this.scatterPlot.camera;
                    const widthHalf = this.props.chartSize.width / 2;
                    const heightHalf = this.props.chartSize.height / 2;
                    const pos = new Vector3();
                    let selectedIndex = -1;
                    const tolerance = 2;
                    if (this.lastHoverIndex !== -1) {
                        pos.x = positions[this.lastHoverIndex * 3];
                        pos.y = positions[this.lastHoverIndex * 3 + 1];
                        pos.z = positions[this.lastHoverIndex * 3 + 2];
                        pos.project(camera);
                        pos.x = (pos.x * widthHalf) + widthHalf;
                        pos.y = -(pos.y * heightHalf) + heightHalf;
                        if (Math.abs(pos.x - point.x) <= tolerance && Math.abs(pos.y - point.y) <= tolerance) {
                            selectedIndex = this.lastHoverIndex;
                        }
                    }

                    if (selectedIndex === -1) {
                        for (let i = 0, j = 0, k = 0; i < trace.npoints; i++, j += 4, k += 3) {
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
                    this.lastHoverIndex = selectedIndex;
                    if (selectedIndex !== -1) {
                        let value = trace.values[selectedIndex];
                        let categoryObject = this.props.categoricalNames[trace.name];
                        if (categoryObject) {
                            let renamedValue = categoryObject[value];
                            if (renamedValue != null) {
                                value = renamedValue;
                            }
                        }

                        if (typeof value === 'number') {
                            value = numberFormat2f(value);
                            if (value.endsWith('.00')) {
                                value = value.substring(0, value.lastIndexOf('.'));
                            }
                        }
                        this.props.setTooltip('' + value);
                    } else {
                        this.props.setTooltip('');
                    }

                }

            };
            this.scatterPlot.lassoCallback = (points, appendToSelection) => {
                const trace = this.props.trace;
                const positions = trace.positions;
                const camera = this.scatterPlot.camera;
                const widthHalf = this.props.chartSize.width / 2;
                const heightHalf = this.props.chartSize.height / 2;
                const pos = new Vector3();
                const selectedIndices = new Set();

                for (let i = 0, j = 0, k = 0; i < trace.npoints; i++, j += 4, k += 3) {
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
                    this.props.onSelected({name: getEmbeddingKey(trace.embedding)});
                } else {
                    this.props.onSelected({
                        name: getEmbeddingKey(trace.embedding),
                        clear: !appendToSelection,
                        value: {basis: trace.embedding, indices: selectedIndices}
                    });
                }
            };
            this.scatterPlot.boxCallback = (rect, appendToSelection) => {
                if (this.scatterPlot.interactionMode === 'PAN') {
                    return;
                }
                const trace = this.props.trace;
                const positions = trace.positions;
                const camera = this.scatterPlot.camera;
                const widthHalf = this.props.chartSize.width / 2;
                const heightHalf = this.props.chartSize.height / 2;
                const pos = new Vector3();
                const selectedIndices = new Set();

                for (let i = 0, j = 0, k = 0; i < trace.npoints; i++, j += 4, k += 3) {
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
                    this.props.onSelected({name: getEmbeddingKey(trace.embedding)});
                } else {
                    this.props.onSelected({
                        name: getEmbeddingKey(trace.embedding),
                        clear: !appendToSelection,
                        value: {basis: trace.embedding, indices: selectedIndices}
                    });
                }
            };
            const canvas = this.containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
            const webglcontextlost = (e) => {
                console.log('lost webgl context');
                e.preventDefault();
                this.scatterPlot = null;
            };
            const webglcontextrestored = (e) => {
                console.log('restored webgl context');
                e.preventDefault();
                this.setState({forceUpdate: !this.state.forceUpdate});
            };
            canvas.addEventListener('webglcontextlost', webglcontextlost);
            canvas.addEventListener('webglcontextrestored', webglcontextrestored);
        }

    }


    draw() {
        const {
            obsCat,
            cachedData,
            trace,
            markerOpacity,
            unselectedMarkerOpacity,
            selection,
            pointSize,
            chartOptions,
            categoricalNames
        } = this.props;
        updateScatterChart(this.scatterPlot, trace, selection, markerOpacity, unselectedMarkerOpacity, pointSize,
            categoricalNames, chartOptions, obsCat, cachedData);
    }


    render() {
        return <>
            <div className={this.props.classes.root}>
                <ChartToolbar
                    dragmode={this.props.chartOptions.dragmode}
                    // editSelection={this.props.chartOptions.editSelection}
                    onGallery={this.props.onGallery}
                    animating={this.props.chartOptions.animating}
                    showFog={this.props.chartOptions.showFog}
                    onShowFog={this.onShowFog}
                    is3d={this.props.trace && this.props.trace.z != null}
                    toggleAnimation={this.onToggleAnimation}
                    onSaveImage={this.onSaveImage}
                    onDragMode={this.onDragMode}
                    onCopyImage={this.onCopyImage}
                    // onEditSelection={this.onEditSelection}
                    onShowAxis={this.onShowAxis}
                    onHome={this.resetCamera}
                    showAxis={this.props.chartOptions.showAxis}
                >
                </ChartToolbar>
            </div>

            <div style={{
                display: 'inline-block',
                width: this.props.chartSize.width,
                height: this.props.chartSize.height
            }}
                 ref={this.containerElementRef}>
            </div>
        </>;
    }
}

export default withStyles(styles)(ScatterChartThree);



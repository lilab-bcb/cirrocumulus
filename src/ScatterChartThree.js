import {Typography} from '@material-ui/core';
import withStyles from '@material-ui/core/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import {bind} from 'lodash';
import React from 'react';
import {Color, Vector3, Vector4} from 'three';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {saveImage} from './ChartUtil';
import {numberFormat} from './formatters';
import {createScatterPlot, getCategoryLabelsPositions, updateScatterChart} from './ThreeUtil';
import {getChartSize, isPointInside} from './util';

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

export const POINT_VISUALIZER_ID = 'SPRITES';
export const LABELS_VISUALIZER_ID = 'CANVAS_LABELS';

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
        this.tooltipElementRef = React.createRef();
        this.scatterPlot = null;
        this.chartSize = getChartSize();
        this.lastHoverIndex = -1;
        // window.addEventListener('resize', () => {
        //     scatterGL.resize();
        // });
    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.draw();
    }

    componentDidMount() {
        this.init();
        this.draw(); // TODO fix 2x draw
        this.draw();
    }

    componentWillUnmount() {
        // if (this.scatterGL != null) {
        //     this.scatterGL.dispose();
        // }
    }


    calculatePointSize(traceInfo) {
        const n = traceInfo.npoints;
        const SCALE = 200;
        const LOG_BASE = 8;
        const DIVISOR = 1.5;
        // Scale point size inverse-logarithmically to the number of points.
        const pointSize = SCALE / Math.log(n) / Math.log(LOG_BASE);
        return traceInfo.dimensions === 3 ? pointSize : pointSize / DIVISOR;
    }

    drawContext(context, chartSize, format) {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, categoricalNames, chartOptions} = this.props;
        const showLabels = this.props.chartOptions.showLabels && traceInfo.isCategorical;

        const pointSize = this.calculatePointSize(traceInfo);
        let scaleFactor = this.props.pointSize;
        const PI2 = 2 * Math.PI;
        const colors = traceInfo.colors;
        const positions = traceInfo.positions;
        const camera = this.scatterPlot.camera;
        const width = chartSize.width;
        const height = chartSize.height;
        const widthHalf = width / 2;
        const heightHalf = height / 2;
        const colorScale = scaleLinear().domain([0, 1]).range([0, 255]);
        const npoints = traceInfo.npoints;

        const is3d = traceInfo.dimensions === 3;
        let outputPointSize = 0;
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
        const pos = new Vector3();
        let cameraSpacePos = new Vector4();
        // vec4 cameraSpacePos = modelViewMatrix * vec4(position, 1.0);
        // Project vertex in camera-space to screen coordinates using the camera's
        // projection matrix.
        // gl_Position = projectionMatrix * cameraSpacePos;


        let object = spriteVisualizer.points;
        let modelViewMatrix = object.modelViewMatrix.clone();
        modelViewMatrix.multiplyMatrices(camera.matrixWorldInverse, object.matrixWorld);
        let gl_PointSize = (outputPointSize * scaleFactor) / 4;
        const showFog = chartOptions.showFog;
        for (let i = 0, j = 0, k = 0; i < npoints; i++, j += 4, k += 3) {
            const isSelected = selection.size === 0 || selection.has(i);
            pos.x = positions[k];
            pos.y = positions[k + 1];
            pos.z = positions[k + 2];
            pos.project(camera);

            let r = colors[j];
            let g = colors[j + 1];
            let b = colors[j + 2];
            let a = isSelected ? markerOpacity : unselectedMarkerOpacity;
            if (is3d) {
                cameraSpacePos.x = positions[k];
                cameraSpacePos.y = positions[k + 1];
                cameraSpacePos.z = positions[k + 2];
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
        if (showLabels) {
            context.fillStyle = 'black';

            let font = format === 'svg' ? 'serif' : 'Roboto Condensed';
            context.font = '24px ' + font; // FIXME adjust size
            const labelsPositions = getCategoryLabelsPositions(traceInfo, categoricalNames);
            for (let i = 0, k = 0; i < labelsPositions.labels.length; i++, k += 3) {
                pos.x = labelsPositions.positions[k];
                pos.y = labelsPositions.positions[k + 1];
                pos.z = labelsPositions.positions[k + 2];
                pos.project(camera);
                pos.x = (pos.x * widthHalf) + widthHalf;
                pos.y = -(pos.y * heightHalf) + heightHalf;
                context.fillText(labelsPositions.labels[i], pos.x, pos.y);
            }
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
        const {traceInfo} = this.props;
        saveImage(traceInfo, this.chartSize, bind(this.drawContext, this), format);
    };

    onEditSelection = () => {
        this.props.chartOptions.editSelection = !this.props.chartOptions.editSelection;
        this.props.setChartOptions(this.props.chartOptions);
    };


    onShowLabels = () => {
        this.props.chartOptions.showLabels = !this.props.chartOptions.showLabels;
        this.props.setChartOptions(this.props.chartOptions);
    };

    onDarkMode = () => {
        this.props.chartOptions.darkMode = !this.props.chartOptions.darkMode;
        this.scatterPlot.scene.background = this.props.chartOptions.darkMode ? new Color("rgb(0, 0, 0)") : null;
        this.props.setChartOptions(this.props.chartOptions);
    };

    onGallery = () => {
        this.props.onGallery();
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
            this.scatterPlot = createScatterPlot(containerElement, false, true);
            const axes = this.scatterPlot.scene.getObjectByName('axes');
            if (axes) {
                axes.visible = this.props.chartOptions.showAxis;
            }
            let spriteVisualizer = getVisualizer(this.scatterPlot, POINT_VISUALIZER_ID);
            spriteVisualizer.styles.fog.enabled = this.props.chartOptions.showFog;
            this.scatterPlot.hoverCallback = (point) => {
                if (point == null) {
                    this.tooltipElementRef.current.innerHTML = ' ';
                } else {
                    const traceInfo = this.props.traceInfo;
                    const positions = traceInfo.positions;
                    const camera = this.scatterPlot.camera;
                    const widthHalf = this.chartSize.width / 2;
                    const heightHalf = this.chartSize.height / 2;
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
                        for (let i = 0, j = 0, k = 0; i < traceInfo.npoints; i++, j += 4, k += 3) {
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
                        let value = traceInfo.values[selectedIndex];
                        let categoryObject = this.props.categoricalNames[traceInfo.name];
                        if (categoryObject) {
                            let renamedValue = categoryObject[value];
                            if (renamedValue != null) {
                                value = renamedValue;
                            }
                        }

                        if (typeof value === 'number') {
                            value = numberFormat(value);
                        }
                        this.tooltipElementRef.current.innerHTML = '' + value;
                    } else {
                        this.tooltipElementRef.current.innerHTML = ' ';
                    }

                }

            };
            this.scatterPlot.lassoCallback = (points) => {
                const traceInfo = this.props.traceInfo;
                const positions = traceInfo.positions;
                const camera = this.scatterPlot.camera;
                const widthHalf = this.chartSize.width / 2;
                const heightHalf = this.chartSize.height / 2;
                const pos = new Vector3();
                const selectedPoints = [];

                for (let i = 0, j = 0, k = 0; i < traceInfo.npoints; i++, j += 4, k += 3) {
                    pos.x = positions[k];
                    pos.y = positions[k + 1];
                    pos.z = positions[k + 2];
                    pos.project(camera);
                    pos.x = (pos.x * widthHalf) + widthHalf;
                    pos.y = -(pos.y * heightHalf) + heightHalf;
                    if (isPointInside(pos, points)) {
                        selectedPoints.push(i);
                    }
                }

                if (selectedPoints.length === 0) {
                    this.props.onDeselect({name: getEmbeddingKey(traceInfo.embedding)});
                } else {
                    this.props.onSelected({
                        name: getEmbeddingKey(traceInfo.embedding),
                        clear: !this.props.chartOptions.editSelection,
                        value: {basis: traceInfo.embedding, points: selectedPoints}
                    });
                }
            };
            this.scatterPlot.boxCallback = (rect) => {
                if (this.scatterPlot.interactionMode === 'PAN') {
                    return;
                }
                const traceInfo = this.props.traceInfo;
                const positions = traceInfo.positions;
                const camera = this.scatterPlot.camera;
                const widthHalf = this.chartSize.width / 2;
                const heightHalf = this.chartSize.height / 2;
                const pos = new Vector3();
                const selectedPoints = [];

                for (let i = 0, j = 0, k = 0; i < traceInfo.npoints; i++, j += 4, k += 3) {
                    pos.x = positions[k];
                    pos.y = positions[k + 1];
                    pos.z = positions[k + 2];
                    pos.project(camera);
                    pos.x = (pos.x * widthHalf) + widthHalf;
                    pos.y = -(pos.y * heightHalf) + heightHalf;
                    if (pos.x >= rect.x && pos.x <= (rect.x + rect.width) && pos.y >= rect.y && pos.y <= (rect.y + rect.height)) {
                        selectedPoints.push(i);
                    }
                }

                if (selectedPoints.length === 0) {
                    this.props.onDeselect({name: getEmbeddingKey(traceInfo.embedding)});
                } else {
                    this.props.onSelected({
                        name: getEmbeddingKey(traceInfo.embedding),
                        clear: !this.props.chartOptions.editSelection,
                        value: {basis: traceInfo.embedding, points: selectedPoints}
                    });
                }
            };
        }
        ;


        const canvas = this.containerElementRef.current.querySelector('canvas');
        canvas.style.outline = '0px';
    }


    draw() {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, pointSize, chartOptions, categoricalNames} = this.props;
        updateScatterChart(this.scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize,
            chartOptions.showLabels, categoricalNames, chartOptions.showFog, chartOptions.showAxis, chartOptions.darkMode);
    }


    render() {
        return <React.Fragment>
            <div className={this.props.classes.root}>
                <ChartToolbar
                    dragmode={this.props.chartOptions.dragmode}
                    animating={this.props.chartOptions.animating}
                    editSelection={this.props.chartOptions.editSelection}
                    showLabels={this.props.chartOptions.showLabels}
                    darkMode={this.props.chartOptions.darkMode}
                    onDarkMode={this.onDarkMode}
                    showFog={this.props.chartOptions.showFog}
                    onShowFog={this.onShowFog}
                    is3d={this.props.traceInfo && this.props.traceInfo.z != null}
                    toggleAnimation={this.onToggleAnimation}
                    onSaveImage={this.onSaveImage}
                    onShowLabels={this.onShowLabels}
                    onDragMode={this.onDragMode}
                    onCopyImage={this.onCopyImage}
                    onEditSelection={this.onEditSelection}
                    onShowAxis={this.onShowAxis}
                    showAxis={this.props.chartOptions.showAxis}
                    onGallery={this.onGallery}>
                </ChartToolbar>
                <Typography color="textPrimary" ref={this.tooltipElementRef} style={{
                    display: 'inline-block',
                    paddingLeft: 5,
                    verticalAlign: 'top'
                }}>&nbsp;</Typography>
            </div>

            <div style={{
                display: 'inline-block',
                width: this.chartSize.width,
                height: this.chartSize.height
            }}
                 ref={this.containerElementRef}>
            </div>
        </React.Fragment>;
    }
}

export default withStyles(styles)(ScatterChartThree);



import withStyles from '@material-ui/core/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import React from 'react';
import {Vector3, Vector4} from 'three';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {drawColorScheme} from './ColorSchemeLegend';
import {numberFormat} from './formatters';
import {drawCategoricalLegend, getCategoricalLegendSize} from './LegendDrawer';
import {createScatterPlot, getCategoryLabelsPositions, updateScatterChart} from './ThreeUtil';
import {getChartSize, isPointInside, setClipboardData} from './util';

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

        this.state = {animating: false, dragmode: 'pan', editSelection: false, showLabels: false};
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
        console.log('unmount');
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
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, categoricalNames} = this.props;
        const showLabels = this.state.showLabels && traceInfo.isCategorical;

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

        let spriteVisualizer;
        for (let i = 0; i < this.scatterPlot.visualizers.length; i++) {
            if (this.scatterPlot.visualizers[i].id === 'SPRITES') {
                spriteVisualizer = this.scatterPlot.visualizers[i];
                break;
            }
        }

        let object = spriteVisualizer.points;
        let modelViewMatrix = object.modelViewMatrix.clone();
        modelViewMatrix.multiplyMatrices(camera.matrixWorldInverse, object.matrixWorld);
        let gl_PointSize = (outputPointSize * scaleFactor) / 4;
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
                const fogDepth = pointSize / outputPointSize * 1.2;
                gl_PointSize = (outputPointSize * scaleFactor) / 4;
                const fogFactor = smoothstep(fog.near, fog.far, fogDepth);
                r = mix(r, fog.color.r, fogFactor);
                g = mix(g, fog.color.g, fogFactor);
                b = mix(b, fog.color.b, fogFactor);
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

    onCopyImage = (event) => {
        const size = this.chartSize;
        const canvas = document.createElement('canvas');
        canvas.width = size.width; // * window.devicePixelRatio;
        canvas.height = size.height; // * window.devicePixelRatio;
        const context = canvas.getContext('2d');
        this.drawContext(context, this.chartSize);
        const url = canvas.toDataURL();
        setClipboardData([{
            format: 'image/png',
            data: '<img src="' + url + '">'
        }]);
    };


    onSaveImage = (format) => {
        const {traceInfo} = this.props;
        let context;
        let canvas = null;
        const chartSize = this.chartSize;
        const totalSize = {width: this.chartSize.width, height: this.chartSize.height};
        let name = traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        if (format !== 'svg') {
            canvas = document.createElement('canvas');
            canvas.width = 100;
            canvas.height = 100;
            context = canvas.getContext('2d');
        } else {
            context = new window.C2S(100, 100);
        }
        if (!traceInfo.continuous) {
            const legendSize = getCategoricalLegendSize(context, name, traceInfo.colorScale.domain());
            totalSize.width += legendSize.width;
            totalSize.height = Math.max(legendSize.height, chartSize.height);
        } else {
            totalSize.height += 150;

        }
        if (format === 'svg') {
            context = new window.C2S(totalSize.width, totalSize.height);
        } else {
            canvas.width = totalSize.width * window.devicePixelRatio;
            canvas.height = totalSize.height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.fillStyle = 'white';
            context.fillRect(0, 0, totalSize.width, totalSize.height);

        }

        this.drawContext(context, chartSize, format);

        if (!traceInfo.continuous) {
            context.translate(chartSize.width, 2);
            drawCategoricalLegend(context, traceInfo.colorScale, name, traceInfo.colorScale.domain());
        } else {
            context.translate(chartSize.width / 2 - 75, chartSize.height + 2);
            drawColorScheme(context, 150, 20, traceInfo.colorScale, true);
        }

        if (format === 'svg') {
            let svg = context.getSerializedSvg();
            // let prefix = [];
            // prefix.push('<?xml version="1.0" encoding="utf-8"?>\n');
            // prefix.push('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"' +
            //     ' "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
            // svg = prefix.join('') + svg;
            let blob = new Blob([svg], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, name + '.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, name + '.png', true);
            });
        }


    };
    onToggleAnimation = () => {
        if (this.scatterPlot.orbitIsAnimating()) {
            this.scatterPlot.stopOrbitAnimation();
            this.setState({animating: false});
        } else {
            this.scatterPlot.startOrbitAnimation();
            this.setState({animating: true});
        }
    };

    onEditSelection = () => {
        this.setState((state, props) => {
            return {editSelection: !state.editSelection};
        });
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
        this.setState({dragmode: mode});
    };

    onShowLabels = () => {
        this.setState((state, props) => {
            return {showLabels: !state.showLabels};
        });

    };

    onGallery = () => {
        this.props.onGallery();
    };

    init() {

        if (this.scatterPlot == null) {
            const containerElement = this.containerElementRef.current;
            this.scatterPlot = createScatterPlot(containerElement, false);
            this.scatterPlot.hoverCallback = (point) => {
                if (point == null) {
                    this.tooltipElementRef.current.innerHTML = ' ';
                } else {
                    const traceInfo = this.props.traceInfo;
                    let value = traceInfo.values[point];
                    if (typeof value === 'number') {
                        value = numberFormat(value);
                    }
                    this.tooltipElementRef.current.innerHTML = '' + value;
                }

            };
            this.scatterPlot.lassoCallback = (points) => {
                const traceInfo = this.props.traceInfo;
                const positions = traceInfo.positions;
                const camera = this.scatterPlot.camera;
                const chartSize = this.chartSize;
                const width = chartSize.width;
                const height = chartSize.height;
                const widthHalf = width / 2;
                const heightHalf = height / 2;
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
                        clear: !this.state.editSelection,
                        value: {basis: traceInfo.embedding, points: selectedPoints}
                    });
                }
            };
            this.scatterPlot.selectCallback = (selectedpoints) => {
                if (this.scatterPlot.interactionMode === 'PAN') {
                    return;
                }
                if (selectedpoints != null && selectedpoints.length === 0) {
                    selectedpoints = null;
                }
                const traceInfo = this.props.traceInfo;
                if (selectedpoints == null) {
                    this.props.onDeselect({name: getEmbeddingKey(traceInfo.embedding)});
                } else {
                    // let xmin = Number.MAX_VALUE;
                    // let ymin = Number.MAX_VALUE;
                    // let zmin = Number.MAX_VALUE;
                    // let xmax = -Number.MAX_VALUE;
                    // let ymax = -Number.MAX_VALUE;
                    // let zmax = -Number.MAX_VALUE;
                    //
                    // const is3d = traceInfo.dimensions === 3;
                    // selectedpoints.forEach(index => {
                    //     const x = traceInfo.x[index];
                    //     xmin = Math.min(xmin, x);
                    //     xmax = Math.max(xmax, x);
                    //     const y = traceInfo.y[index];
                    //     ymin = Math.min(ymin, y);
                    //     ymax = Math.max(ymax, y);
                    //     if (is3d) {
                    //         const z = traceInfo.z[index];
                    //         zmin = Math.min(zmin, z);
                    //         zmax = Math.max(zmax, z);
                    //     }
                    // });
                    //
                    //
                    // let path = {
                    //     shape: 'rect',
                    //     x: xmin,
                    //     y: ymin,
                    //     width: xmax - xmin,
                    //     height: ymax - ymin
                    // };
                    // if (is3d) {
                    //     path.shape = 'rect 3d';
                    //     if (zmin > zmax) {
                    //         const tmp = zmin;
                    //         zmin = zmax;
                    //         zmax = tmp;
                    //     }
                    //     path.z = zmin;
                    //     path.depth = zmax - zmin;
                    // }
                    //
                    // this.props.onSelected({
                    //     name: getEmbeddingKey(traceInfo.embedding),
                    //     clear: !this.state.editSelection,
                    //     value: {basis: traceInfo.embedding, path: path}
                    // });
                    this.props.onSelected({
                        name: getEmbeddingKey(traceInfo.embedding),
                        clear: !this.state.editSelection,
                        value: {basis: traceInfo.embedding, points: selectedpoints}
                    });
                }
            };


            const canvas = this.containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
        }
    }

    draw() {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, color, pointSize, categoricalNames} = this.props;
        updateScatterChart(this.scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize, this.state.showLabels, categoricalNames);
    }


    render() {
        return <React.Fragment>
            <div className={this.props.classes.root}>
                <ChartToolbar
                    dragmode={this.state.dragmode}
                    animating={this.state.animating}
                    editSelection={this.state.editSelection}
                    showLabels={this.state.showLabels}
                    is3d={this.props.traceInfo && this.props.traceInfo.z != null}
                    onHome={this.onHome}
                    toggleAnimation={this.onToggleAnimation}
                    onSaveImage={this.onSaveImage}
                    onShowLabels={this.onShowLabels}
                    onDragMode={this.onDragMode}
                    onCopyImage={this.onCopyImage}
                    onEditSelection={this.onEditSelection}
                    onGallery={this.onGallery}>
                </ChartToolbar>
                <div ref={this.tooltipElementRef} style={{
                    display: 'inline-block',
                    paddingLeft: 5,
                    verticalAlign: 'top'
                }}>&nbsp;</div>
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



import {Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import PhotoLibraryIcon from '@material-ui/icons/PhotoLibrary';
import {scaleLinear} from 'd3-scale';
import React from 'react';
import {Vector3, Vector4} from 'three';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {numberFormat} from './formatters';
import {createScatterPlot} from './ThreeUtil';
import {getChartSize, setClipboardData} from './util';

export function updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize) {
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
    scatterPlot.setPointColors(colors);
    scatterPlot.setPointPositions(positions);
    scatterPlot.setDimensions(traceInfo.dimensions);
    // const {scaleDefault, scaleSelected, scaleHover} = this.scatterPlot.styles.point;

    const scale = new Float32Array(traceInfo.npoints);
    scale.fill(pointSize);
    scatterPlot.setPointScaleFactors(scale);
    scatterPlot.render();
}

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

class ScatterChartThree extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();
        this.tooltipElementRef = React.createRef();
        this.scatterPlot = null;
        this.chartSize = getChartSize();

        this.state = {animating: false, dragmode: 'pan', editSelection: false};
        // window.addEventListener('resize', () => {
        //     scatterGL.resize();
        // });

    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.draw();

    }

    componentDidMount() {
        this.init();
        this.draw();
        this.draw();
    }

    componentWillUnmount() {
        console.log('unmount');
        // if (this.scatterGL != null) {
        //     this.scatterGL.dispose();
        // }
    }

    onCopyImage = (event) => {
        const canvas = this.containerElementRef.current.querySelector('canvas');
        const url = canvas.toDataURL();
        setClipboardData([{
            format: 'text/html',
            data: '<img src="' + url + '">'
        }], true);

    };


    calculatePointSize(traceInfo) {
        const n = traceInfo.npoints;
        const SCALE = 200;
        const LOG_BASE = 8;
        const DIVISOR = 1.5;
        // Scale point size inverse-logarithmically to the number of points.
        const pointSize = SCALE / Math.log(n) / Math.log(LOG_BASE);
        return traceInfo.dimensions === 3 ? pointSize : pointSize / DIVISOR;
    }

    drawContext(context, chartSize) {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection} = this.props;
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

        let object = this.scatterPlot.visualizers.get('SPRITES').points;
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
    }

    savePng() {
        const {traceInfo} = this.props;
        const size = this.chartSize;
        const canvas = document.createElement('canvas');
        canvas.width = size.width; // * window.devicePixelRatio;
        canvas.height = size.height; // * window.devicePixelRatio;
        const context = canvas.getContext('2d');
        this.drawContext(context, this.chartSize);
        // context.scale(window.devicePixelRatio, window.devicePixelRatio);
        let name = traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        canvas.toBlob(blob => {
            window.saveAs(blob, name + '.png', true);
        });


    }

    saveSvg() {
        const {traceInfo} = this.props;
        let context = new window.C2S(this.chartSize.width, this.chartSize.height);
        this.drawContext(context, this.chartSize);
        let svg = context.getSerializedSvg();
        let blob = new Blob([svg], {
            type: 'text/plain;charset=utf-8'
        });
        let name = traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        window.saveAs(blob, name + '.svg');
    }

    onSaveImage = (format) => {
        // if (this.scatterPlot.orbitIsAnimating()) {
        //     this.scatterPlot.stopOrbitAnimation();
        //     this.setState({animating: false});
        // }
        if (format === 'svg') {
            this.saveSvg();
        } else {
            this.savePng();
        }

        // const renderer = new SVGRenderer();
        // document.body.appendChild(renderer.domElement);
        // renderer.setSize(window.innerWidth, window.innerHeight);
        // renderer.setQuality('low');
        // renderer.render(this.scatterGL.scatterPlot.scene, this.scatterGL.scatterPlot.camera);
        // console.log(renderer.domElement);

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
        }
        this.setState({dragmode: mode});
    };

    onGallery = (event) => {
        event.preventDefault();
        this.props.onGallery();
    };

    init() {

        if (this.scatterPlot == null) {
            const containerElement = this.containerElementRef.current;
            this.scatterPlot = createScatterPlot(containerElement);
            this.scatterPlot.hoverCallback = (point) => {
                if (point == null) {
                    this.tooltipElementRef.current.innerHTML = ' ';
                } else {
                    const traceInfo = this.props.traceInfo;
                    let value = traceInfo.values[point];
                    if (typeof value === 'number') {
                        value = numberFormat(value);
                    }
                    this.tooltipElementRef.current.innerHTML = ' ' + value;
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

                    let xmin = Number.MAX_VALUE;
                    let ymin = Number.MAX_VALUE;
                    let zmin = Number.MAX_VALUE;
                    let xmax = -Number.MAX_VALUE;
                    let ymax = -Number.MAX_VALUE;
                    let zmax = -Number.MAX_VALUE;
                    const is3d = traceInfo.dimensions === 3;
                    selectedpoints.forEach(index => {
                        const x = traceInfo.x[index];
                        xmin = Math.min(xmin, x);
                        xmax = Math.max(xmax, x);
                        const y = traceInfo.y[index];
                        ymin = Math.min(ymin, y);
                        ymax = Math.max(ymax, y);
                        if (is3d) {
                            const z = traceInfo.z[index];
                            zmin = Math.min(zmin, z);
                            zmax = Math.max(zmax, z);
                        }
                    });


                    let path = {shape: 'rect', x: xmin, y: ymin, width: xmax - xmin, height: ymax - ymin};
                    if (is3d) {
                        path.shape = 'rect 3d';
                        path.z = zmin;
                        path.depth = zmax - zmin;
                    }
                    this.props.onSelected({
                        name: getEmbeddingKey(traceInfo.embedding),
                        clear: !this.state.editSelection,
                        value: {basis: traceInfo.embedding, path: path}
                    });
                }
            };

            // this.scatterGL = new ScatterGL(this.containerElementRef.current, {
            //     renderMode: 'POINT',
            //     rotateOnStart: false,
            //     showLabelsOnHover: false,
            //
            //     onHover: (point) => {
            //
            //     }
            // });


            const canvas = this.containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
        }
    }

    draw() {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, color, pointSize} = this.props;
        updateScatterChart(this.scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize);
    }


    render() {
        return <React.Fragment>
            <ChartToolbar
                dragmode={this.state.dragmode}
                animating={this.state.animating}
                editSelection={this.state.editSelection}
                is3d={this.props.traceInfo && this.props.traceInfo.z != null}
                onHome={this.onHome}
                toggleAnimation={this.onToggleAnimation}
                onSaveImage={this.onSaveImage}
                onDragMode={this.onDragMode}
                onCopyImage={this.onCopyImage}
                onEditSelection={this.onEditSelection}>
            </ChartToolbar>
            <Tooltip title={"View Gallery"}>
                <IconButton edge={false} size={'small'}

                            aria-label="View Gallery" onClick={this.onGallery}>
                    <PhotoLibraryIcon/>
                </IconButton>
            </Tooltip>


            <div ref={this.tooltipElementRef} style={{
                display: 'inline-block',
                textOverflow: 'hidden',
                paddingLeft: 5
            }}>&nbsp;
            </div>

            <div style={{
                display: 'inline-block',
                position: 'relative',
                width: this.chartSize.width,
                height: this.chartSize.height
            }}
                 ref={this.containerElementRef}>

            </div>
        </React.Fragment>;
    }
}

export default ScatterChartThree;




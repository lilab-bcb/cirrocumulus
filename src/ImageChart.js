import CircularProgress from '@mui/material/CircularProgress';
import withStyles from '@mui/styles/withStyles';
import {bind, uniqueId} from 'lodash';
import OpenSeadragon from 'openseadragon';
import React from 'react';
import simplify from 'simplify-js';
import {getEmbeddingKey} from './actions';
import CanvasOverlayHd from './CanvasOverlayHd';
import ChartToolbar from './ChartToolbar';
import {saveImage} from './ChartUtil';
import {numberFormat2f} from './formatters';
import OpenseadragonSvgOverlay from './OpenseadragonSvgOverlay';
import {getCategoryLabelsPositions, getLabels} from './ThreeUtil';
import {arrayToSvgPath, isPointInside} from './util';
import {DEFAULT_STYLE, SCATTER_TRANSITION, setTooltipPosition} from './CirroTooltip';


export function getSpotRadius(trace, pointSize) {
    return pointSize * (trace.embedding.spatial.spot_diameter ? trace.embedding.spatial.spot_diameter / 2 : 20);
}

export function drawEmbeddingImage(context, chartSize, trace, selection, markerOpacity, unselectedMarkerOpacity, chartOptions, categoricalNames, obsCat, cachedData, spotRadius) {
    if (trace.tileSource.ready) {
        const img = trace.tileSource.levels[trace.tileSource.levels.length - 1].context2D.canvas;
        if (chartSize == null) {
            chartSize = {width: img.width, height: img.height};
        }
        const zoom = Math.min(chartSize.width / img.width, chartSize.height / img.height);
        context.drawImage(img, 0, 0, img.width * zoom, img.height * zoom);
        context.scale(zoom, zoom);
        drawSpots(context, zoom, trace, selection, markerOpacity, unselectedMarkerOpacity, spotRadius);
        drawLabels(context, zoom, trace, chartOptions, categoricalNames, obsCat, cachedData);
        context.setTransform(1, 0, 0, 1, 0, 0);
    }
}

function drawLabels(context, zoom, trace, chartOptions, categoricalNames, obsCat, cachedData) {
    const showLabels = obsCat.length > 0;
    if (showLabels) {
        context.textAlign = 'center';
        context.textBaseline = "middle";
        const darkMode = true; // assume image is dark chartOptions.darkMode;
        const fontSize = Math.ceil(chartOptions.labelFontSize / zoom);
        context.fillStyle = darkMode ? 'white' : 'black';
        context.strokeStyle = darkMode ? 'rgba(0,0,0,0.9)' : 'rgba(255,255,255,0.9)';
        context.lineWidth = chartOptions.labelStrokeWidth;

        context.font = fontSize + 'px Roboto Condensed,Helvetica,Arial,sans-serif';

        const labelsPositions = getCategoryLabelsPositions(trace.embedding, obsCat, cachedData);
        const labels = getLabels(obsCat, labelsPositions.labels, categoricalNames);
        for (let i = 0, index = 0, n = labels.length; i < n; i++, index += 3) {
            let x = labelsPositions.positions[index];
            let y = labelsPositions.positions[index + 1];
            context.strokeText('' + labels[i], x, y);
            context.fillText('' + labels[i], x, y);
        }
    }
}

function drawSpots(context, zoom, trace, selection, markerOpacity, unselectedMarkerOpacity, spotRadius) {
    context.lineWidth = 2 * 1 / zoom;
    if (context.setLineDash) {
        context.setLineDash([2, 2]);
    }
    const isSelectionEmpty = selection == null;
    const indices = trace.indices;
    if (!isSelectionEmpty) { // draw unselected cells 1st
        context.globalAlpha = unselectedMarkerOpacity;
        for (let i = 0; i < trace.x.length; i++) {
            let index = indices[i];
            let x = trace.x[index];
            let y = trace.y[index];
            if (!selection.has(index)) {
                context.fillStyle = trace.colors[index];
                context.beginPath();
                context.arc(x, y, spotRadius, 0, Math.PI * 2, true);
                context.fill();
            }
        }
        context.globalAlpha = markerOpacity;
        for (let i = 0; i < trace.x.length; i++) {
            let index = indices[i];
            let x = trace.x[index];
            let y = trace.y[index];
            if (selection.has(index)) {
                context.fillStyle = trace.colors[index];
                context.beginPath();
                context.arc(x, y, spotRadius, 0, Math.PI * 2, true);
                context.fill();
            }
        }
    } else {
        context.globalAlpha = markerOpacity;
        for (let i = 0; i < trace.x.length; i++) {
            let index = indices[i];
            let x = trace.x[index];
            let y = trace.y[index];
            context.fillStyle = trace.colors[index];
            context.beginPath();
            context.arc(x, y, spotRadius, 0, Math.PI * 2, true);
            context.fill();
        }
    }
    context.globalAlpha = 1;
    if (context.setLineDash) {
        context.setLineDash([]);
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
        position: 'absolute',
        zIndex: 1,
        top: 0,
        left: 0,
        display: 'inline-block',
        verticalAlign: 'top',
        whiteSpace: 'nowrap',
        overflow: 'hidden'
    }
});

class ImageChart extends React.PureComponent {

    constructor(props) {
        super(props);
        this.id = uniqueId('cirro-image');
        this.editSelection = false;
        this.state = {loading: false};
        this.tooltipElement = null;
    }

    findPointsInPolygon(points) {
        let data = this.props.trace;
        // let spotRadius = data[0].embedding.spotDiameter / 2;
        const indices = [];
        for (let i = 0; i < data.x.length; i++) {
            if (isPointInside({x: data.x[i], y: data.y[i]}, points)) {
                indices.push(i);
            }
        }
        return indices;
    }

    findPointsInRectangle(rect) {
        const data = this.props.trace;
        const spotRadius = getSpotRadius(data, this.props.pointSize);
        const indices = [];
        const x = parseFloat(rect.getAttribute('x'));
        const y = parseFloat(rect.getAttribute('y'));
        const x2 = x + parseFloat(rect.getAttribute('width'));
        const y2 = y + parseFloat(rect.getAttribute('height'));
        for (let i = 0; i < data.x.length; i++) {
            const px = data.x[i];
            const py = data.y[i];
            if (px <= x2 && x <= px + spotRadius && py <= y2 && y <= py + spotRadius) {
                indices.push(i);
            }
        }
        return indices;
    }


    findPointIndex(xpix, ypix) {
        const data = this.props.trace;
        const spotRadius = getSpotRadius(data, this.props.pointSize);
        // x is spot center
        for (let i = 0; i < data.x.length; i++) {
            if (Math.abs(data.x[i] - xpix) <= spotRadius && Math.abs(data.y[i] - ypix) <= spotRadius) {
                return i;
            }
        }
        return -1;
    }

    setTooltip(point, webPoint) {
        if (point == null) {
            this.tooltipElement.innerHTML = '';
            setTooltipPosition(this.tooltipElement, -1, -1, SCATTER_TRANSITION);
        } else {
            const trace = this.props.trace;
            let value = trace.values[point];
            const categoryObject = this.props.categoricalNames[trace.name];
            if (categoryObject) {
                let renamedValue = categoryObject[value];
                if (renamedValue != null && renamedValue.newValue != null) {
                    value = renamedValue.newValue;
                }
            }
            if (typeof value === 'number') {
                value = numberFormat2f(value);
                if (value.endsWith('.00')) {
                    value = value.substring(0, value.lastIndexOf('.'));
                }
            }
            const parentRect = this.tooltipElement.parentElement.getBoundingClientRect();
            this.tooltipElement.innerHTML = '' + value;
            const tooltipRect = this.tooltipElement.getBoundingClientRect();
            let left = webPoint.x + 8;
            if ((left + tooltipRect.width + 4) >= parentRect.width) {
                left = webPoint.x - 8 - tooltipRect.width;
            }
            let top = webPoint.y + 8;
            if ((top + tooltipRect.height + 4) >= parentRect.height) {
                top = webPoint.y - 8 - tooltipRect.height;
            }

            setTooltipPosition(this.tooltipElement, left + 'px', top + 'px', SCATTER_TRANSITION);
        }
    }


    drawContext(context, chartSize) {
        const img = this.viewer.source.levels[this.viewer.source.levels.length - 1].context2D.canvas;
        if (chartSize == null) {
            chartSize = {width: img.width, height: img.height};
        }

        const zoom = Math.min(chartSize.width / img.width, chartSize.height / img.height);
        context.drawImage(img, 0, 0, img.width * zoom, img.height * zoom);
        this._drawOverlay({context: context, zoom: zoom});
    }

    _drawOverlay(opts) {
        const context = opts.context;
        const trace = this.props.trace;
        const selection = this.props.selection;
        const markerOpacity = this.props.markerOpacity;
        const unselectedMarkerOpacity = this.props.unselectedMarkerOpacity;
        const spotRadius = getSpotRadius(trace, this.props.pointSize);
        drawSpots(context, opts.zoom, trace, selection, markerOpacity, unselectedMarkerOpacity, spotRadius);
        drawLabels(context, opts.zoom, trace, this.props.chartOptions, this.props.categoricalNames, this.props.obsCat, this.props.cachedData);
    }

    createViewer() {
        if (this.viewer) {
            this.viewer.destroy();
        }
        // let tileSource = new OpenSeadragon.ImageTileSource({
        //     url: url,
        //     buildPyramid: true,
        //     crossOriginPolicy: "Anonymous"
        // });
        if (!this.props.trace.tileSource.ready) {
            this.setState({loading: true});
            this.props.trace.tileSource.addOnceHandler('ready', (src) => {
                this.setState({loading: false});
            });
        } else {
            this.setState({loading: false});
        }
        this.viewer = new OpenSeadragon({
            id: this.id,
            gestureSettingsMouse: {dblClickToZoom: true, clickToZoom: false},
            autoResize: true,
            showFullPageControl: false,
            collectionMode: false,
            // visibilityRatio: 0.2,
            showNavigationControl: false,
            // prefixUrl: 'https://cdn.jsdelivr.net/npm/openseadragon@2.4/build/openseadragon/images/',
            tileSources: this.props.trace.tileSource
        });
        const viewer = this.viewer;
        const _this = this;
        this.canvasOverlay = new CanvasOverlayHd(this.viewer, {
            onRedraw: function (opts) {
                _this._drawOverlay(opts);
            }
        });
        const tooltipElement = document.createElement("div");
        for (const key in DEFAULT_STYLE) {
            tooltipElement.style[key] = DEFAULT_STYLE[key];
        }
        this.viewer.canvas.appendChild(tooltipElement);
        this.tooltipElement = tooltipElement;

        let lassoPathArray = [];
        let startCoordinates = [0, 0];
        let lastBoundingBox = {x: 0, y: 0, width: 0, height: 0};
        const tracker = new OpenSeadragon.MouseTracker({
            element: this.viewer.container,
            moveHandler: function (event) {
                if (viewer.world == null) {
                    _this.setTooltip();
                } else if (_this.props.chartOptions.dragmode === 'pan' && viewer.world.getItemCount() > 0) {
                    const webPoint = event.position;
                    const viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                    const imagePoint = viewer.viewport.viewportToImageCoordinates(viewportPoint);
                    const point = _this.findPointIndex(imagePoint.x, imagePoint.y);
                    if (point !== -1) {
                        _this.setTooltip(point, webPoint);
                    } else {
                        _this.setTooltip();
                    }
                } else {
                    _this.setTooltip();
                }
            }
        });
        tracker.setTracking(true);

        this.viewer.addHandler('canvas-double-click', function (event) {
            if (_this.props.chartOptions.dragmode === 'pan') {
                const webPoint = event.position;
                const viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                const imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
                const point = _this.findPointIndex(imagePoint.x, imagePoint.y);
                if (point !== -1) {
                    event.preventDefaultAction = true;
                    _this.props.handleClick({
                        name: _this.props.trace.name,
                        value: _this.props.trace.values[point],
                        shiftKey: false,
                        metaKey: false
                    });
                }
            }
        });


        this.viewer.addHandler('canvas-drag', function (event) {
            if ((_this.props.chartOptions.dragmode === 'lasso' || _this.props.chartOptions.dragmode === 'select') && viewer.world.getItemCount() > 0) {
                event.preventDefaultAction = true;
                const webPoint = event.position;
                const viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                const imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
                lassoPathArray.push({x: imagePoint.x, y: imagePoint.y});
                if (_this.props.chartOptions.dragmode === 'lasso') {
                    lassoPathArray = simplify(lassoPathArray);
                    lassoPath.setAttribute('d', arrayToSvgPath(lassoPathArray));
                } else {
                    lastBoundingBox.x = Math.min(imagePoint.x, startCoordinates[0]);
                    lastBoundingBox.y = Math.max(imagePoint.y, startCoordinates[1]);
                    lastBoundingBox.width =
                        Math.max(imagePoint.x, startCoordinates[0]) - lastBoundingBox.x;
                    lastBoundingBox.height =
                        lastBoundingBox.y - Math.min(imagePoint.y, startCoordinates[1]);

                    rectElement.setAttribute('x', '' + lastBoundingBox.x);
                    rectElement.setAttribute(
                        'y',
                        '' + (lastBoundingBox.y - lastBoundingBox.height)
                    );
                    rectElement.setAttribute('width', '' + lastBoundingBox.width);
                    rectElement.setAttribute('height', '' + lastBoundingBox.height);
                }
            }
        });

        this.viewer.addHandler('canvas-press', function (event) {
            if (_this.props.chartOptions.dragmode === 'lasso' || _this.props.chartOptions.dragmode === 'select') {
                event.preventDefaultAction = true;
                _this.editSelection = event.metaKey || event.ctrlKey;
                const webPoint = event.position;
                const viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                const imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);

                if (_this.props.chartOptions.dragmode === 'lasso') {
                    lassoPathArray = [];
                    lassoPathArray.push({x: imagePoint.x, y: imagePoint.y});
                    lassoPathArray = simplify(lassoPathArray);
                    lassoPath.setAttribute('d', arrayToSvgPath(lassoPathArray));
                } else {
                    startCoordinates = [imagePoint.x, imagePoint.y];
                    lastBoundingBox = {
                        x: startCoordinates[0],
                        y: startCoordinates[1],
                        width: 1,
                        height: 1
                    };

                }
            }
        });

        this.viewer.addHandler('canvas-release', function (event) {
            if (_this.props.chartOptions.dragmode === 'lasso') {
                event.preventDefaultAction = true;
                const indices = new Set(_this.findPointsInPolygon(lassoPathArray));
                lassoPathArray = [];
                lassoPath.setAttribute('d', '');
                _this.props.onSelected({
                    name: getEmbeddingKey(_this.props.trace.embedding),
                    clear: !_this.editSelection,
                    value: {basis: _this.props.trace.embedding, indices: indices}
                });
            } else if (_this.props.chartOptions.dragmode === 'select') {
                event.preventDefaultAction = true;
                const indices = new Set(_this.findPointsInRectangle(rectElement));
                rectElement.removeAttribute('x');
                rectElement.removeAttribute('y');
                rectElement.removeAttribute('width');
                rectElement.removeAttribute('height');
                _this.props.onSelected({
                    name: getEmbeddingKey(_this.props.trace.embedding),
                    clear: !_this.editSelection,
                    value: {basis: _this.props.trace.embedding, indices: indices}
                });
            }
        });

        viewer.addHandler('canvas-exit', function (event) {
            _this.setTooltip();
        });

        const svgOverlay = new OpenseadragonSvgOverlay(viewer);
        const lassoPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        const rectElement = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        [lassoPath, rectElement].forEach(elem => {
            elem.setAttribute('stroke', 'black');
            elem.setAttribute('stroke-width', '3px');
            elem.setAttribute('fill', '#0bb');
            elem.setAttribute('fill-opacity', '0.3');
            elem.setAttribute('stroke-dasharray', '2 2');
            svgOverlay.node().appendChild(elem);
        });
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (this.canvasOverlay) {
            this.canvasOverlay._updateCanvas();
        }
        // if (this.viewer && prevProps.chartSize !== this.props.chartSize) {
        //     this.viewer.viewport.resize();
        // }
    }

    componentWillUnmount() {
        if (this.viewer) {
            this.viewer.destroy();
        }
        this.viewer = null;
    }

    componentDidMount() {
        if (this.viewer == null) {
            this.createViewer();
        }
    }


    onSaveImage = (format) => {
        const {trace} = this.props;
        const img = this.viewer.source.levels[this.viewer.source.levels.length - 1].context2D.canvas;
        saveImage(trace, {width: img.width, height: img.height}, bind(this.drawContext, this), format);
    };


    onZoomIn = () => {
        this.viewer.viewport.zoomBy(this.viewer.zoomPerClick / 1.0);
        this.viewer.viewport.applyConstraints();
    };

    onZoomOut = () => {
        this.viewer.viewport.zoomBy(1.0 / this.viewer.zoomPerClick);
        this.viewer.viewport.applyConstraints();
    };

    onHome = () => {
        this.viewer.viewport.goHome();
        this.viewer.viewport.applyConstraints();
    };

    onDragMode = (mode) => {
        this.props.chartOptions.dragmode = mode;
        this.props.setChartOptions(this.props.chartOptions);
    };


    render() {
        return <>
            <div className={this.props.classes.root}>
                <ChartToolbar
                    dragmode={this.props.chartOptions.dragmode}
                    // editSelection={this.props.chartOptions.editSelection}
                    onGallery={this.props.onGallery}
                    animating={false}
                    onZoomIn={this.onZoomIn}
                    onZoomOut={this.onZoomOut}
                    is3d={false}
                    onHome={this.onHome}
                    onSaveImage={this.onSaveImage}
                    onDragMode={this.onDragMode}
                    // onEditSelection={this.onEditSelection}
                >
                </ChartToolbar>
            </div>
            <div style={{
                display: 'inline-block',
                width: this.props.chartSize.width,
                height: this.props.chartSize.height
            }}
                 id={this.id}>
            </div>
            {this.state.loading && <CircularProgress style={{
                position: 'absolute',
                left: this.props.chartSize.width / 2,
                top: this.props.chartSize.height / 2
            }} size={20}/>}
        </>;
    }
}


export default withStyles(styles)(ImageChart);

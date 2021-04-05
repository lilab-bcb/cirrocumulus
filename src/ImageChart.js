import CircularProgress from '@material-ui/core/CircularProgress';
import withStyles from '@material-ui/core/styles/withStyles';
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
        const darkMode = true;
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
    const isSelectionEmpty = selection.size === 0;
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
};
const styles = theme => ({

    root: {
        '& > *': {
            margin: theme.spacing(.4),
        },
        '& > .MuiIconButton-root': {
            padding: 0,
        },
        position: 'absolute',
        zIndex: 1,
        top: 0,
        left: 0,
        display: 'inline-block',
        verticalAlign: 'top',
        whiteSpace: 'nowrap',
        overflow: 'hidden',
    }
});

class ImageChart extends React.PureComponent {

    constructor(props) {
        super(props);
        this.id = uniqueId('cirro-image');
        this.editSelection = false;
        this.state = {loading: false};
    }

    findPointsInPolygon(points) {
        let data = this.props.trace;
        // let spotRadius = data[0].embedding.spotDiameter / 2;
        let indices = [];
        for (let i = 0; i < data.x.length; i++) {
            if (isPointInside({x: data.x[i], y: data.y[i]}, points)) {
                indices.push(i);
            }
        }
        return indices;
    }

    findPointsInRectangle(rect) {
        let data = this.props.trace;
        const spotRadius = getSpotRadius(data, this.props.pointSize);
        let indices = [];
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
        let data = this.props.trace;
        const spotRadius = getSpotRadius(data, this.props.pointSize);
        for (let i = 0; i < data.x.length; i++) {
            if (Math.abs(data.x[i] - xpix) <= spotRadius && Math.abs(data.y[i] - ypix) <= spotRadius) {
                return i;
            }
        }
        return -1;
    }

    setTooltip(xpix, ypix) {
        let trace = this.props.trace;
        const point = this.findPointIndex(xpix, ypix);
        if (point != -1) {
            let value = trace.values[point];
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
        let context = opts.context;
        let trace = this.props.trace;
        const selection = this.props.selection;
        let markerOpacity = this.props.markerOpacity;
        let unselectedMarkerOpacity = this.props.unselectedMarkerOpacity;
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
        let viewer = this.viewer;

        let _this = this;
        this.canvasOverlay = new CanvasOverlayHd(this.viewer, {
            onRedraw: function (opts) {
                _this._drawOverlay(opts);
            },
        });
        // let tooltip = document.createElement("div");
        // tooltip.style.background = 'rgba(0,0,0,0.5)';
        // tooltip.style.color = 'white';
        // tooltip.style.position = 'absolute';
        // this.viewer.canvas.appendChild(tooltip);

        let lassoPathArray = [];
        let startCoordinates = [0, 0];
        let lastBoundingBox = {x: 0, y: 0, width: 0, height: 0};

        this.viewer.innerTracker.moveHandler = function (event) {
            if (viewer.world == null) {
                _this.props.setTooltip('');
            } else if (_this.props.chartOptions.dragmode === 'pan' && viewer.world.getItemCount() > 0) {
                let webPoint = event.position;
                let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
                _this.setTooltip(imagePoint.x, imagePoint.y);
                // if (title !== '') {
                //     tooltip.innerHTML = title;
                //     tooltip.style.left = (webPoint.x + 8) + 'px';
                //     tooltip.style.top = (webPoint.y + 8) + 'px';
                //     tooltip.style.display = '';
                // } else {
                //     tooltip.style.display = 'none';
                // }
            } else {
                _this.props.setTooltip('');
            }
        };
        //
        // this.viewer.innerTracker.scrollHandler = function (event) {
        //     if (lassoState > -1) {
        //         event.preventDefaultAction = true;
        //     }
        // };


        this.viewer.addHandler('canvas-drag', function (event) {
            if ((_this.props.chartOptions.dragmode === 'lasso' || _this.props.chartOptions.dragmode === 'select') && viewer.world.getItemCount() > 0) {
                event.preventDefaultAction = true;
                let webPoint = event.position;
                let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
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
                let webPoint = event.position;
                let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);

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
                        height: 1,
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

        // this.viewer.innerTracker.clickHandler = function (event) {
        //     if (_this.props.chartOptions.dragmode === 'pan') {
        //         let webPoint = event.position;
        //         let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
        //         let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
        //         const point = _this.findPointIndex(imagePoint.x, imagePoint.y);
        //         if (point === -1) {
        //             //   this.props.onSelected({name: getEmbeddingKey(trace.embedding)});
        //         } else {
        //             _this.props.onSelected({
        //                 name: getEmbeddingKey(_this.props.trace.embedding),
        //                 clear: !_this.props.chartOptions.editSelection,
        //                 value: {basis: _this.props.trace.embedding, points: [point]}
        //             });
        //         }
        //     }
        // };

        viewer.addHandler('canvas-exit', function (event) {
            _this.props.setTooltip('');

        });

        let svgOverlay = new OpenseadragonSvgOverlay(viewer);
        let lassoPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        let rectElement = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
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
            this.createViewer(this.props.trace.url);
        }
    }


    onSaveImage = (format) => {
        const {trace} = this.props;
        const img = this.viewer.source.levels[this.viewer.source.levels.length - 1].context2D.canvas;
        saveImage(trace, {width: img.width, height: img.height}, bind(this.drawContext, this), format);
    };

    // onEditSelection = () => {
    //     this.props.chartOptions.editSelection = !this.props.chartOptions.editSelection;
    //     this.props.setChartOptions(this.props.chartOptions);
    // };


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

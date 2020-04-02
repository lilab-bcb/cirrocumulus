import {uniqueId} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import CanvasOverlayHd from './CanvasOverlayHd';
import {intFormat} from './formatters';
import OpenseadragonSvgOverlay from './OpenseadragonSvgOverlay';
import {arrayToSvgPath, isPointInside} from './util';

const OpenSeadragon = window.OpenSeadragon;


class ImageChart extends React.PureComponent {

    constructor(props) {
        super(props);
        this.id = uniqueId('cirro-image');
    }

    getSpotsInPolygon(points) {
        let data = this.props.data;
        // let spotRadius = data[0].embedding.spotDiameter / 2;
        let indices = [];
        for (let i = 0; i < data[0].x.length; i++) {
            if (isPointInside({x: data[0].x[i], y: data[0].y[i]}, points)) {
                indices.push(i);
            }
        }
        return indices;
    }

    spotMouseMove(xpix, ypix) {
        let data = this.props.data;
        let spotRadius = data[0].embedding.spotDiameter / 2;
        let coordIndex = -1;
        for (let i = 0; i < data[0].x.length; i++) {
            if (Math.abs(data[0].x[i] - xpix) <= spotRadius && Math.abs(data[0].y[i] - ypix) <= spotRadius) {
                coordIndex = i;
                break;
            }
        }
        let title = [];
        if (coordIndex != -1) {
            title.push(data[0].text[coordIndex]);
            title.push('x: ' + intFormat(Math.round(data[0].x[coordIndex])));
            title.push('y: ' + intFormat(Math.round(data[0].y[coordIndex])));
        }
        return title.join('<br/>');
    }


    drawSpots = (opts) => {
        let context = opts.context;
        let data = this.props.data;
        let spotRadius = data[0].embedding.spotDiameter / 2;
        let selectedOpacity = data[0].marker.opacity;
        let unselectedOpacity = data[0].unselected.marker.opacity;

        let selectedPoints = data[0].selectedpoints;
        let selectedPointsIndex = 0;
        context.globalAlpha = selectedOpacity;
        context.lineWidth = 2 * 1 / opts.zoom;
        context.setLineDash([2, 2]);

        // all spots
        if (data[0].name === '__count' && !this.props.data.continuous) {
            context.strokeStyle = 'black';
            for (let i = 0; i < data[0].x.length; i++) {
                let x = data[0].x[i];
                let y = data[0].y[i];
                let isSelected = false;
                if (selectedPoints == null || i === selectedPoints[selectedPointsIndex]) {
                    isSelected = true;
                    selectedPointsIndex++;
                }
                context.globalAlpha = isSelected ? selectedOpacity : unselectedOpacity;
                context.beginPath();
                context.arc(x, y, spotRadius, 0, Math.PI * 2, true);
                context.stroke();
            }
        } else {
            for (let i = 0; i < data[0].x.length; i++) {
                let x = data[0].x[i];
                let y = data[0].y[i];
                let isSelected = false;
                if (selectedPoints == null || i === selectedPoints[selectedPointsIndex]) {
                    isSelected = true;
                    selectedPointsIndex++;
                }
                context.globalAlpha = isSelected ? selectedOpacity : unselectedOpacity;
                context.fillStyle = data[0].colors[i];
                context.beginPath();
                context.arc(x, y, spotRadius, 0, Math.PI * 2, true);
                context.fill();
            }
        }
        context.globalAlpha = 1;
        context.setLineDash([]);
    };

    createViewer(url) {
        let viewer = new OpenSeadragon({
            id: this.id,
            gestureSettingsMouse: {dblClickToZoom: true, clickToZoom: false},
            autoResize: false,
            showFullPageControl: false,
            collectionMode: false,
            // visibilityRatio: 0.2,
            prefixUrl: 'https://cdn.jsdelivr.net/npm/openseadragon@2.4/build/openseadragon/images/',
            tileSources: {
                type: 'image',
                url: url,
                buildPyramid: true,
            },
        });

        this.viewer = viewer;
        let _this = this;


        this.canvasOverlay = new CanvasOverlayHd(this.viewer, {
            onRedraw: function (opts) {
                _this.drawSpots(opts);
            },
        });
        let tooltip = document.createElement("div");
        tooltip.style.background = 'rgba(0,0,0,0.5)';
        tooltip.style.color = 'white';
        tooltip.style.position = 'absolute';
        this.viewer.canvas.appendChild(tooltip);

        let selectionMode = -1; // -1=no lasso, 0=start, 1=dragging
        let lassoPathArray = [];


        this.viewer.innerTracker.moveHandler = function (event) {
            if (selectionMode === -1 && viewer.world.getItemCount() > 0) {
                let webPoint = event.position;
                let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
                let title = _this.spotMouseMove(imagePoint.x, imagePoint.y);
                if (title !== '') {
                    tooltip.innerHTML = title;
                    tooltip.style.left = (webPoint.x + 8) + 'px';
                    tooltip.style.top = (webPoint.y + 8) + 'px';
                    tooltip.style.display = '';
                } else {
                    tooltip.style.display = 'none';
                }
            } else {
                tooltip.style.display = 'none';
            }
        };
        //
        // this.viewer.innerTracker.scrollHandler = function (event) {
        //     if (selectionMode > -1) {
        //         event.preventDefaultAction = true;
        //     }
        // };


        this.viewer.addHandler('canvas-drag', function (event) {
            if (selectionMode > -1 && viewer.world.getItemCount() > 0) {
                event.preventDefaultAction = true;
                if (selectionMode === 0) {
                    lassoPathArray = [];
                }
                let webPoint = event.position;
                let viewportPoint = viewer.viewport.pointFromPixel(webPoint);
                let imagePoint = viewer.world.getItemAt(0).viewportToImageCoordinates(viewportPoint, true);
                lassoPathArray.push({x: imagePoint.x, y: imagePoint.y});
                lassoPath.setAttribute('d', arrayToSvgPath(lassoPathArray));
                selectionMode = 1;
            }
        });

        this.viewer.innerTracker.clickHandler = function (event) {
            if (selectionMode > -1) {
                selectionMode = -1;
                event.preventDefaultAction = true;
                let selectedpoints = _this.getSpotsInPolygon(lassoPathArray);
                if (selectedpoints.length === 0) {
                    selectedpoints = null;
                }
                _this.props.data[0].selectedpoints = selectedpoints;
                if (selectedpoints == null) {
                    _this.props.onDeselect();
                } else {
                    _this.props.onSelected({points: [{data: _this.props.data[0]}]});
                }
                lassoPathArray = [];
                lassoPath.setAttribute('d', '');
            }

        };

        viewer.addHandler('canvas-exit', function (event) {
            tooltip.style.display = 'none';
        });

        let svgOverlay = new OpenseadragonSvgOverlay(viewer);
        let lassoPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        lassoPath.setAttribute('stroke', 'black');
        lassoPath.setAttribute('stroke-width', '2px');
        lassoPath.setAttribute('fill', '#0bb');
        lassoPath.setAttribute('fill-opacity', '0.1');
        lassoPath.setAttribute('stroke-dasharray', '2 2');

        svgOverlay.node().appendChild(lassoPath);
        let pathOverlay = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (this.canvasOverlay) {
            this.canvasOverlay._updateCanvas();
        }
    }

    componentWillUnmount() {
    }

    componentDidMount() {
        if (this.viewer == null) {
            if (this.props.data[0].url == null) {
                this.props.data[0].image.then(response => {
                    return response.blob();
                }).then(blob => URL.createObjectURL(blob)).then(url => {
                    this.props.data[0].url = url;
                    this.createViewer(url);
                });
            } else {
                this.createViewer(this.props.data[0].url);
            }
        }
    }

    render() {
        return <div style={{display: 'inline-block', width: this.props.layout.width, height: this.props.layout.height}}
                    id={this.id}/>;
    }
}


const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        embeddingChartSize: state.embeddingChartSize
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(ImageChart));
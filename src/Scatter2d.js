import {scaleLinear} from 'd3-scale';
import {clientPoint} from 'd3-selection';
import {saveAs} from 'file-saver';
import React from 'react';
import ChartToolbar from './ChartToolbar';
import {arrayToSvgPath, isPointInside} from './PlotUtil';

class Scatter2d extends React.PureComponent {

    constructor(props) {
        super(props);
        this.chartRef = React.createRef();
        this.tooltipRef = React.createRef();
        this.lassoRef = React.createRef();
        this.brushRef = React.createRef();
        this.backingScale = 1;
        this.interactionMode = ChartToolbar.MODE_LASSO;
        this.isPointerDown = false;
        this.selectionMode = -1; // -1=no lasso/brush, 0=start, 1=dragging
        this.lassoPathArray = [];
        this.translateX = 0;
        this.translateY = 0;
    }

    componentWillUnmount() {
        window.removeEventListener('mouseup', this.onMouseUp);
    }

    onHome = () => {
        this.translateX = 0;
        this.translateY = 0;
        this.chartRef.current.style.transform = 'translate(' + this.translateX + 'px,' + this.translateY + 'px)';
    };

    onSaveImage = () => {
        let layout = this.props.layout;
        let height = layout.height;
        let width = layout.width;
        let context = new window.C2S(width, height);
        this.drawContext(context);
        let svg = context.getSerializedSvg();
        let blob = new Blob([svg], {
            type: 'text/plain;charset=utf-8'
        });
        saveAs(blob, 'image.svg');
    };

    onLasso = () => {
        this.interactionMode = ChartToolbar.MODE_LASSO;
    };

    onZoomIn = () => {
        this.interactionMode = ChartToolbar.MODE_ZOOM_IN;
    };

    onZoomOut = () => {
        this.interactionMode = ChartToolbar.MODE_ZOOM_OUT;
    };

    onPan = () => {
        this.interactionMode = ChartToolbar.MODE_PAN;
    };

    onBrush = () => {
        this.interactionMode = ChartToolbar.MODE_BRUSH;
    };

    onMouseUp = (event) => {
        this.isPointerDown = false;
        window.removeEventListener('mouseup', this.onMouseUp);

        if (this.interactionMode === ChartToolbar.MODE_LASSO) {
            let trace = this.props.data[0];

            let selectedpoints = [];

            let point = {x: 0, y: 0};

            for (let i = 0, n = trace.x.length; i < n; i++) {
                let x = this.xToPix(trace.x[i]); // includes translate
                let y = this.yToPix(trace.y[i]);
                point.x = x;
                point.y = y;
                if (isPointInside(point, this.lassoPathArray)) {
                    selectedpoints.push(i);
                }
            }
            if (selectedpoints.length === 0) {
                selectedpoints = null;
            }
            trace.selectedpoints = selectedpoints;
            if (selectedpoints == null) {
                this.props.onDeselect();
            } else {
                this.props.onSelected({points: [{data: this.props.data[0]}]});
            }
            this.lassoPathArray = [];
            this.lassoRef.current.setAttribute('d', '');
        }

    };

    onMouseDown = (event) => {

        if (event.button === 0) {
            window.addEventListener('mouseup', this.onMouseUp);
            this.isPointerDown = true;
            if (this.interactionMode === ChartToolbar.MODE_LASSO) {
                this.lassoPathArray = [];
                let coords = clientPoint(event.target, event);
                this.lassoPathArray.push({x: coords[0], y: coords[1]});
                this.lassoRef.current.setAttribute('d', arrayToSvgPath(this.lassoPathArray));
            } else if (this.interactionMode === ChartToolbar.MODE_PAN) {
                this.mouseDownPoint = {x: event.clientX - this.translateX, y: event.clientY - this.translateY};
            }
        }
    };

    onMouseMove = (event) => {

        if (this.isPointerDown) {
            if (this.interactionMode === ChartToolbar.MODE_LASSO) {
                let coords = clientPoint(event.target, event);
                this.lassoPathArray.push({x: coords[0], y: coords[1]});
                this.lassoRef.current.setAttribute('d', arrayToSvgPath(this.lassoPathArray));
            } else if (this.interactionMode === ChartToolbar.MODE_PAN) {
                this.translateX = event.clientX - this.mouseDownPoint.x;
                this.translateY = event.clientY - this.mouseDownPoint.y;
                this.chartRef.current.style.transform = 'translate(' + this.translateX + 'px,' + this.translateY + 'px)';
            }
        }
    };


    redraw() {

        //let quadTree = d3.geom.quadtree(data);
        let backingScale = this.backingScale;
        let node = this.chartRef.current;
        const context = node.getContext('2d');
        let layout = this.props.layout;
        let height = layout.height;
        let width = layout.width;
        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);
        this.drawContext(context);
    }


    xToPix(value) {
        return this.xToPixScale(value) + this.translateX;
    }

    yToPix(value) {
        return this.yToPixScale(value) + this.translateY;
    }

    pixToX(pix) {
        return this.xToPixScale.invert()(pix - this.translateX);
    }

    pixToY(pix) {
        return this.yToPixScale.invert()(pix - this.translateY);
    }


    drawContext(context) {
        let data = this.props.data;
        let layout = this.props.layout;
        let height = layout.height;
        let width = layout.width;
        let trace = data[0];
        if (trace.xmin == null) {
            let xmin = Number.MAX_VALUE;
            let xmax = -Number.MAX_VALUE;
            let ymin = Number.MAX_VALUE;
            let ymax = -Number.MAX_VALUE;
            for (let i = 0, n = trace.x.length; i < n; i++) {
                let x = trace.x[i];
                let y = trace.y[i];
                xmin = x < xmin ? x : xmin;
                xmax = x > xmax ? x : xmax;
                ymin = y < ymin ? y : ymin;
                ymax = y > ymax ? y : ymax;
            }
            trace.xmin = xmin;
            trace.xmax = xmax;
            trace.ymin = ymin;
            trace.ymax = ymax;
        }
        let markerSize = trace.marker.size;
        let unselectedMarkerSize = trace.unselected.marker.size;
        let maxMarkerSize = Math.max(markerSize, unselectedMarkerSize);
        let markerOpacity = trace.marker.opacity;
        let xToPixScale = scaleLinear().domain([trace.xmin, trace.xmax]).range([maxMarkerSize, width - maxMarkerSize]);
        let yToPixScale = scaleLinear().domain([trace.ymin, trace.ymax]).range([height - maxMarkerSize, maxMarkerSize]);
        this.xToPixScale = xToPixScale;
        this.yToPixScale = yToPixScale;
        const PI2 = 2 * Math.PI;
        if (trace.selectedpoints == null) {
            context.globalAlpha = markerOpacity;
            for (let i = 0, n = trace.x.length; i < n; i++) {
                let xpix = this.xToPixScale(trace.x[i]);
                let ypix = this.yToPixScale(trace.y[i]);
                context.fillStyle = trace.marker.color[i];
                context.beginPath();
                context.arc(xpix, ypix, markerSize, 0, PI2);
                context.closePath();
                context.fill();
            }
        } else {
            let unselectedOpacity = trace.unselected.marker.opacity;

            let selectedPoints = trace.selectedpoints;
            let selectedPointsIndex = 0;
            for (let i = 0, n = trace.x.length; i < n; i++) {
                let xpix = this.xToPixScale(trace.x[i]);
                let ypix = this.yToPixScale(trace.y[i]);
                let isSelected = false;
                if (i === selectedPoints[selectedPointsIndex]) {
                    isSelected = true;
                    selectedPointsIndex++;
                }
                context.globalAlpha = isSelected ? markerOpacity : unselectedOpacity;
                context.fillStyle = trace.marker.color[i];
                context.beginPath();
                context.arc(xpix, ypix, isSelected ? markerSize : unselectedMarkerSize, 0, PI2);
                context.closePath();
                context.fill();
            }
        }
        context.setTransform(1, 0, 0, 1, 0, 0);

    }

    componentDidMount() {
        this.redraw();
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.data != this.props.data || prevProps.layout != this.props.layout) {
            this.redraw();
        }
    }

    render() {

        let backingScale = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            backingScale = window.devicePixelRatio;
        }
        this.backingScale = backingScale;
        let layout = this.props.layout;
        let height = layout.height;
        let width = layout.width;

        let style = {width: width, height: height, position: 'absolute'};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }

        return (

            <div style={{position: 'relative', width: width, height: height + 16, overflow: 'hidden'}}>

                <canvas width={width * backingScale} height={height * backingScale}
                        ref={this.chartRef}
                        style={style}

                ></canvas>
                <div style={{
                    color: 'white',
                    position: 'absolute',
                    background: 'rgba(0,0,0,0.5)'
                }}
                     ref={this.tooltipRef}></div>
                <svg style={{width: width, height: height, position: 'absolute'}}
                     onMouseMove={this.onMouseMove}
                     onMouseDown={this.onMouseDown}
                >
                    <g>
                        <path fillOpacity={0.1} fill='#0bb' strokeDasharray='2 2' strokeWidth='2px'
                              stroke='black'
                              ref={this.lassoRef}></path>
                        <rect fillOpacity={0.1} fill='#0bb' strokeDasharray='2 2' strokeWidth='2px'
                              stroke='black'
                              ref={this.brushRef}></rect>
                    </g>
                </svg>
                <ChartToolbar onHome={this.onHome}
                              onZoomIn={this.onZoomIn}
                              onZoomOut={this.onZoomOut}
                              onPan={this.onPan}
                              onBrush={this.onBrush}
                              onLasso={this.onLasso}
                              onSaveImage={this.onSaveImage}
                />
            </div>

        );

    }
}

export default Scatter2d;



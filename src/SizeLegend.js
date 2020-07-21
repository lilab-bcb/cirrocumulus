import {scaleLinear} from 'd3-scale';
import React from 'react';
import {numberFormat0} from './formatters';

export function drawSizeLegend(context, scale, nsteps, width, margin = 20) {
    let domain = scale.domain();
    let value = domain[0];

    let stepSize = (domain[1] - domain[0]) / nsteps;
    let legendHeight = 20;

    let valueToX = scaleLinear().range([margin, width - margin]).domain([0, nsteps - 1]).clamp(true);
    let valueToRadius = scaleLinear().range([2, 10]).domain(domain).clamp(true);


    context.textBaseline = 'top';
    context.fillStyle = 'black';
    context.textAlign = 'center';

    for (let i = 0; i < nsteps; i++) {
        if (i === (nsteps - 1)) {
            value = domain[1];
        }
        let pix = valueToX(i);

        let radius = valueToRadius(value);
        context.beginPath();
        context.arc(pix, 10, radius, 0, Math.PI * 2);
        context.stroke();

        context.fillText(numberFormat0(100 * value), pix, legendHeight + 2);

        value += stepSize;
    }
}

class SizeLegend extends React.PureComponent {

    constructor(props) {
        super(props);
        this.ref = React.createRef();
        this.backingScale = 1;
    }

    redraw() {
        let backingScale = this.backingScale;
        let node = this.ref.current;
        const context = node.getContext('2d');
        const height = this.props.height;
        const width = this.props.width;
        context.font = '12px Roboto Condensed,Helvetica,Arial,sans-serif';
        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);
        const scale = this.props.scale;
        if (scale == null) {
            return;
        }

        drawSizeLegend(context, scale, this.props.nsteps || 3, width);
        context.setTransform(1, 0, 0, 1, 0, 0);

    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    componentDidMount() {
        this.redraw();
    }

    render() {

        let backingScale = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            backingScale = window.devicePixelRatio;
        }
        this.backingScale = backingScale;
        let height = this.props.height;
        let width = this.props.width;
        let style = {width: width, height: height};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }
        return (
            <canvas width={width * backingScale} height={height * backingScale} ref={this.ref} style={style}></canvas>);

    }
}

export default SizeLegend;


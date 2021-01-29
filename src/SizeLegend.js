import {scaleLinear} from 'd3-scale';
import React from 'react';
import {CANVAS_FONT} from './ChartUtil';
import {numberFormat2f} from './formatters';
import {stripTrailingZeros} from './util';

export function drawSizeLegend(context, scale, nsteps, width, margin = 20, textColor = 'black') {
    let domain = scale.domain();

    let stepSize = (domain[1] - domain[0]) / (nsteps + 1);
    let steps = [];
    steps.push(domain[0]);
    let value = domain[0];
    for (let i = 0; i < nsteps; i++) {
        value += stepSize;
        steps.push(value);
    }
    steps.push(domain[1]);
    let legendHeight = 20;
    let valueToX = scaleLinear().range([margin, width - margin]).domain([0, steps.length - 1]).clamp(true);
    let radiusDomain = [1, 9];
    if (scale.range()[0] > scale.range()[1]) {
        radiusDomain = [9, 1];
    }
    let valueToRadius = scaleLinear().range(radiusDomain).domain(domain).clamp(true);

    context.textBaseline = 'top';
    context.fillStyle = textColor;
    context.strokeStyle = textColor;
    context.textAlign = 'center';

    for (let i = 0; i < steps.length; i++) {
        let pix = valueToX(i);
        let radius = valueToRadius(steps[i]);
        context.beginPath();
        context.arc(pix, 10, radius, 0, Math.PI * 2);
        context.stroke();
        let text = stripTrailingZeros(numberFormat2f(steps[i]));
        context.fillText(text, pix, legendHeight + 2);

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
        const textColor = this.props.textColor || 'black';
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);
        const scale = this.props.scale;
        if (scale == null) {
            return;
        }

        drawSizeLegend(context, scale, this.props.nsteps || 3, width, 20, textColor);
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


import {scaleLinear} from 'd3-scale';
import React from 'react';
import {numberFormat} from './formatters';

export function drawColorScheme(context, colorScale, textColor = 'black', label = true, width = 150, height = 12) {
    let domain = colorScale.domain();
    if (domain[0] === domain[1]) {
        return;
    }

    let legendHeight = height;
    let gradient = context.createLinearGradient(0, 0, width, legendHeight);
    let valueToFraction = scaleLinear().range([0, 1]).domain(domain).clamp(true);
    gradient.addColorStop(valueToFraction(domain[0]), colorScale(domain[0]));
    gradient.addColorStop(valueToFraction(domain[1]), colorScale(domain[1]));
    context.fillStyle = gradient;
    context.fillRect(0, 0, width, legendHeight);
    context.strokeStyle = 'LightGrey';
    context.strokeRect(0, 0, width, legendHeight);
    if (label) {

        context.textBaseline = 'top';
        context.fillStyle = textColor;

        context.textAlign = 'left';
        context.fillText(numberFormat(domain[0]), 0, legendHeight + 2);
        context.textAlign = 'right';
        context.fillText(numberFormat(domain[1]), width, legendHeight + 2);
    }
}

class ColorSchemeLegend extends React.PureComponent {

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
        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);

        const colorScale = this.props.scale;
        if (colorScale == null) {
            return;
        }

        context.font = '12px Roboto Condensed,Helvetica,Arial,sans-serif';
        drawColorScheme(context, colorScale, textColor, this.props.label, width, height);
        context.setTransform(1, 0, 0, 1, 0, 0);
    }

    componentDidMount() {
        this.redraw();
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    render() {
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            this.backingScale = window.devicePixelRatio;
        }

        let height = this.props.height;
        let width = this.props.width;
        let style = {width: width, height: height};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }
        return (
            <canvas width={width * this.backingScale} height={height * this.backingScale} ref={this.ref}
                    style={style}></canvas>
        );

    }
}

export default ColorSchemeLegend;


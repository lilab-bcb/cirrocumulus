import {scaleLinear} from 'd3-scale';
import React from 'react';
import {numberFormat} from './formatters';

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

        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);
        const colorScale = this.props.scale;
        if (colorScale == null) {
            return;
        }
        let domain = colorScale.domain();
        if (domain[0] === domain[1]) {
            return;
        }
        let value = domain[0];
        let nsteps = this.props.nsteps || 10;

        let stepSize = (domain[1] - domain[0]) / nsteps;
        let legendHeight = 20;
        let gradient = context.createLinearGradient(0, 0, width, legendHeight);
        let valueToFraction = scaleLinear().range([0, 1]).domain(domain).clamp(true);
        for (let i = 0; i < nsteps; i++) {
            if (i === (nsteps - 1)) {
                value = domain[1];
            }
            let f = valueToFraction(value);
            if (!isNaN(f)) {
                let color = colorScale(value);
                gradient.addColorStop(f, color);
            }
            value += stepSize;
        }
        context.fillStyle = gradient;
        context.fillRect(0, 0, width, legendHeight);
        context.strokeStyle = 'LightGrey';
        context.strokeRect(0, 0, width, legendHeight);
        if (this.props.label > 0) {

            context.font = '12px Roboto Condensed,Helvetica,Arial,sans-serif';
            context.textBaseline = 'top';
            context.fillStyle = 'black';

            context.textAlign = 'left';
            context.fillText(numberFormat(domain[0]), 0, legendHeight + 2);
            context.textAlign = 'right';
            context.fillText(numberFormat(domain[1]), width, legendHeight + 2);
        }
        context.setTransform(1, 0, 0, 1, 0, 0);

    }

    componentDidMount() {
        this.redraw();
    }

    handleStart = (e) => {

    };

    handleDrag = (e, data) => {


    };

    handleStop = (e) => {

    };

    componentDidUpdate(prevProps, prevState, snapshot) {
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
        let style = {width: width, height: height, display: 'block'};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }

        return (
            <canvas width={width * backingScale} height={height * backingScale} ref={this.ref} style={style}></canvas>

        );

    }
}

export default ColorSchemeLegend;


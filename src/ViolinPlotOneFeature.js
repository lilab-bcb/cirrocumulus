import Box from '@material-ui/core/Box';
import Typography from '@material-ui/core/Typography';
import {scaleLinear} from 'd3-scale';
import {throttle} from 'lodash';
import React from 'react';
import {CANVAS_FONT} from './ChartUtil';
import {CHIP_SIZE} from './DotPlotCanvas';

export function drawFeature(context, size, feature, data, colorScale, options, drawCategories, categoryColorScales, textColor) {
    const {violinScale, violinHeight, violinWidth, violinShowBoxplot} = options;
    const names = data.map(array => array[0].name); // array of arrays
    const features = data[0].map(item => item.feature);
    const featureIndex = features.indexOf(feature);
    let xmin = Number.MAX_VALUE;
    let xmax = -Number.MAX_VALUE;
    let ymax = -Number.MAX_VALUE;
    for (let i = 0; i < names.length; i++) {
        const item = data[i][featureIndex];
        xmin = Math.min(xmin, item.density.x[0]);
        xmax = Math.max(xmax, item.density.x[item.density.x.length - 1]);
        ymax = Math.max(ymax, item.density.max);
    }

    let yscale;
    if (violinScale === 'area') {
        yscale = scaleLinear().domain([-ymax, ymax]).range([4, violinWidth - 4]); // horizontal position
    }
    const xscale = scaleLinear().domain([xmin, xmax]).range([violinHeight - 1, 1]).nice(); // vertical position
    context.strokeStyle = textColor;
    const boxplotWidth = 6;
    const lineCap = context.lineCap;
    const lineJoin = context.lineJoin;
    for (let categoryIndex = 0; categoryIndex < names.length; categoryIndex++) {
        context.save();
        const item = data[categoryIndex][featureIndex];
        if (violinScale === 'width') {
            yscale = scaleLinear().domain([-item.density.max, item.density.max]).range([4, violinWidth - 4]); // horizontal position
        }
        context.lineCap = 'round';
        context.lineJoin = 'round';
        const density = item.density;
        // context.fillStyle = colorScale(item.mean);
        context.translate(size.x + categoryIndex * violinWidth, 0);
        context.beginPath();
        context.moveTo(yscale(density.y[0]), xscale(density.x[0]));

        for (let i = 1, n = density.x.length; i < n; i++) {
            context.lineTo(yscale(-density.y[i]), xscale(density.x[i]));
        }
        for (let i = density.x.length - 2; i > 0; i--) {
            context.lineTo(yscale(density.y[i]), xscale(density.x[i]));
        }
        // context.closePath();

        // context.fill();
        context.stroke();
        context.lineCap = lineCap;
        context.lineJoin = lineJoin;
        if (violinShowBoxplot) {
            // iqr box
            context.strokeRect(violinWidth / 2 - boxplotWidth / 2, xscale(item.boxplotStats.q3), boxplotWidth, xscale(item.boxplotStats.q1) - xscale(item.boxplotStats.q3));

            // median
            context.beginPath();
            context.moveTo(violinWidth / 2 - boxplotWidth / 2, xscale(item.boxplotStats.median));
            context.lineTo(violinWidth / 2 - boxplotWidth / 2 + boxplotWidth, xscale(item.boxplotStats.median));
            context.stroke();

            // mean
            // context.setLineDash([2, 5]);
            // context.beginPath();
            // context.moveTo(violinWidth / 2 - boxplotWidth / 2, xscale(item.boxplotStats.mean));
            // context.lineTo(violinWidth / 2 - boxplotWidth / 2 + boxplotWidth, xscale(item.boxplotStats.mean));
            // context.stroke();
            // context.setLineDash([]);

            // line from q3 to upperAdjacentValue
            context.beginPath();
            context.moveTo(violinWidth / 2, xscale(item.boxplotStats.upperAdjacentValue));
            context.lineTo(violinWidth / 2, xscale(item.boxplotStats.q3));
            context.stroke();
            // line from q1 to lowerAdjacentValue
            context.beginPath();
            context.moveTo(violinWidth / 2, xscale(item.boxplotStats.q1));
            context.lineTo(violinWidth / 2, xscale(item.boxplotStats.lowerAdjacentValue));
            context.stroke();
        }
        context.restore();
    }


    context.textAlign = "right";
    context.textBaseline = "middle";
    context.fillStyle = textColor;
    context.strokeStyle = textColor;

    const tickWidth = 4;
    let textWidth = size.x - tickWidth;
    // y axis
    context.lineWidth = 0.5;
    context.beginPath();
    context.moveTo(textWidth + tickWidth, xscale(xscale.domain()[0]));
    context.lineTo(textWidth + tickWidth, xscale(xscale.domain()[1]));
    context.stroke();
    // x axis
    context.beginPath();
    context.moveTo(textWidth + tickWidth, xscale(xscale.domain()[0]));
    context.lineTo(size.width + size.x - 4, xscale(xscale.domain()[0]));
    context.stroke();

    context.lineWidth = 1;
    const format = xscale.tickFormat(4);
    const ticks = xscale.ticks(4);
    ticks.forEach(tick => {
        const pix = xscale(tick);
        context.fillText(format(tick), textWidth, pix);
        context.beginPath();
        context.moveTo(textWidth, pix);
        context.lineTo(textWidth + tickWidth, pix);
        context.stroke();
    });

    if (drawCategories) {
        context.textBaseline = 'middle';


        context.textAlign = "left";
        const height = violinHeight + size.y - 1;
        for (let i = 0; i < names.length; i++) {
            const item = data[i][featureIndex];
            const pix = i * violinWidth + violinWidth / 2;
            const name = names[i];
            for (let j = 0; j < name.length; j++) {
                // chip, 2px, text, 4px, ...
                const chipStartCoord = j === 0 ? 0 : size.endCoordinates[j - 1];
                const categoryColorScale = categoryColorScales[j];
                const category = item.categories[j];
                context.fillStyle = categoryColorScale(category);
                context.beginPath();
                context.rect(size.x + pix - CHIP_SIZE + 4, height - CHIP_SIZE - chipStartCoord, CHIP_SIZE, CHIP_SIZE);
                context.fill();
                context.stroke();

                context.save();
                context.fillStyle = textColor;
                context.translate(size.x + pix, height - chipStartCoord - CHIP_SIZE - 2);
                context.rotate(-Math.PI / 2);
                context.fillText(name[j], 0, 0);
                context.restore();
            }
        }
    }

}

export default class ViolinPlotOneFeature extends React.PureComponent {
    constructor(props) {
        super(props);
        this.initialized = false;
        this.canvasRef = React.createRef();
        this.mousemove = throttle(this.mousemove, 500);
    }


    mousemove = (event) => {
        const node = event.target;
        const rect = node.getBoundingClientRect();
        const {data, feature, size, options} = this.props;
        const {violinHeight, violinWidth} = options;
        let xy = [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
        if (xy[0] < size.x || xy[1] >= violinHeight) {
            this.props.onTooltip(null);
        } else {
            const col = Math.floor((xy[0] - size.x) / violinWidth);
            const categories = data.map(array => array[0].name);

            if (col >= 0 && col < categories.length) {
                const features = data[0].map(item => item.feature);
                const featureIndex = features.indexOf(feature);
                this.props.onTooltip(data[col][featureIndex]);
            } else {
                this.props.onTooltip(null);
            }
        }
    };

    mouseout = (event) => {
        this.mousemove.cancel();
        this.props.onTooltip(null);
    };

    redraw() {
        const {categoryColorScales, colorScale, data, options, feature, size, textColor} = this.props;
        if (data == null) {
            return;
        }
        const {violinHeight} = options;
        let devicePixelRatio = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            devicePixelRatio = window.devicePixelRatio;
        }

        const canvas = this.canvasRef.current;
        if (!this.initialized) {
            canvas.addEventListener("mousemove", this.mousemove);
            canvas.addEventListener("mouseout", this.mouseout);
            this.initialized = true;
        }

        const height = violinHeight + size.y;
        const width = size.width + size.x;
        const context = canvas.getContext('2d');
        canvas.width = width * devicePixelRatio;
        canvas.height = height * devicePixelRatio;
        canvas.style.width = width + 'px';
        canvas.style.height = height + 'px';
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        drawFeature(context, size, feature, data, colorScale, options, true, categoryColorScales, textColor);

    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    componentDidMount() {
        this.redraw();
    }

    render() {
        const feature = this.props.feature;
        return (
            <Box borderColor="text.primary" border={1}
                 style={{display: 'inline-block', margin: 2}}>
                <Typography color="textPrimary" component={"h4"}
                            style={{
                                marginTop: '3.2px',
                            }}>{feature}</Typography>
                <canvas ref={this.canvasRef}></canvas>
            </Box>);

    }

}





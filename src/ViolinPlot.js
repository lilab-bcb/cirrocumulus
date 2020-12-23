import {Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';
import PhotoCameraIcon from '@material-ui/icons/PhotoCamera';
import {scaleLinear} from 'd3-scale';
import React from 'react';
import {drawColorScheme} from './ColorSchemeLegend';

const canvasFont = '12px Roboto Condensed,Helvetica,Arial,sans-serif';
const svgFont = '12px Helvetica,Arial,sans-serif';
export default class ViolinPlot extends React.PureComponent {

    constructor(props) {
        super(props);
        this.divRef = React.createRef();
        this.tooltipElementRef = React.createRef();
        this.canvas = null;
        this.state = {saveImageEl: null};
    }


    handleSaveImageMenu = (event) => {
        this.setState({saveImageEl: event.currentTarget});
    };
    handleSaveImageMenuClose = (event) => {
        this.setState({saveImageEl: null});
    };
    handleSaveImage = (format) => {
        this.setState({saveImageEl: null});
        let context;

        let canvas;
        if (format === 'svg') {
            context = new window.C2S(10, 10);
            context.font = svgFont;
        } else {
            canvas = document.createElement('canvas');
            context = canvas.getContext('2d');
            context.font = canvasFont;
        }

        const size = this.getSize(context);

        const colorScaleHeight = 15;
        const height = size.height + size.y + colorScaleHeight + 20;
        const width = Math.max(200, size.width + size.x);
        let scale = 1;

        if (format === 'svg') {
            context = new window.C2S(width, height);
            context.font = svgFont;
        } else {
            canvas.width = width * window.devicePixelRatio;
            canvas.height = height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.font = canvasFont;
        }
        const textColor = 'black';
        // const textColor = this.props.textColor;
        context.fillStyle = textColor === 'white' ? 'black' : 'white';
        context.fillRect(0, 0, width, height);
        this.drawContext(context, size);
        context.translate(4, (size.height + size.y + 4));
        drawColorScheme(context, this.props.colorScale, textColor);

        if (format === 'svg') {
            let svg = context.getSerializedSvg();
            let blob = new Blob([svg], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, this.props.data[0][0].dimension + '.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, this.props.data[0][0].dimension + '.png', true);
            });
        }
    };

    redraw() {
        if (this.props.data == null) {
            return <div/>;
        }
        let devicePixelRatio = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            devicePixelRatio = window.devicePixelRatio;
        }

        if (this.canvas == null) {

            this.canvas = document.createElement('canvas');
            // this.canvas.addEventListener("mousemove", onMouseMove);
            // this.canvas.addEventListener("mouseout", onMouseOut);
            this.divRef.current.append(this.canvas);
        }

        const height = this.size.height + this.size.y;
        const width = this.size.width + this.size.x;
        let canvas = this.canvas;
        const context = canvas.getContext('2d');
        canvas.width = width * devicePixelRatio;
        canvas.height = height * devicePixelRatio;
        canvas.style.width = width + 'px';
        canvas.style.height = height + 'px';
        context.font = canvasFont;

        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        this.drawContext(context, this.size);
    }

    update() {
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        context.font = canvasFont;
        this.size = this.getSize(context);
    }

    drawContext(context, size) {
        const {data} = this.props;
        const features = data[0].map(item => item.feature);
        for (let i = 0; i < features.length; i++) {
            this.drawFeature(context, size, features[i]);
        }
        context.setTransform(1, 0, 0, 1, 0, 0);
    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    componentDidMount() {
        this.redraw();
    }

    drawFeature(context, size, feature) {
        const {data, colorScale, options} = this.props;
        const {violinScale, violinHeight, violinWidth} = options;
        const categories = data.map(array => array[0].name);
        const features = data[0].map(item => item.feature);
        const featureIndex = features.indexOf(feature);
        let xmin = Number.MAX_VALUE;
        let xmax = -Number.MAX_VALUE;
        let ymax = -Number.MAX_VALUE;
        for (let i = 0; i < categories.length; i++) {
            const item = data[i][featureIndex];
            xmin = Math.min(xmin, item.density.x[0]);
            xmax = Math.max(xmax, item.density.x[item.density.x.length - 1]);
            ymax = Math.max(ymax, item.density.max);
        }


        let yscale;
        if (violinScale === 'area') {
            yscale = scaleLinear().domain([-ymax, ymax]).range([4, violinWidth - 4]); // horizontal position
        }
        const xscale = scaleLinear().domain([xmin, xmax]).range([violinHeight - 10, 10]).nice(); // vertical position
        context.strokeStyle = 'black';

        for (let categoryIndex = 0; categoryIndex < categories.length; categoryIndex++) {
            context.save();
            const item = data[categoryIndex][featureIndex];
            if (violinScale === 'width') {
                yscale = scaleLinear().domain([-item.density.max, item.density.max]).range([4, violinWidth - 4]); // horizontal position
            }

            const density = item.density;
            context.fillStyle = colorScale(item.mean);
            context.translate(size.x + categoryIndex * violinWidth, violinHeight * featureIndex);
            context.beginPath();
            context.moveTo(yscale(density.y[0]), xscale(density.x[0]));

            for (let i = 1, n = density.x.length; i < n; i++) {
                context.lineTo(yscale(-density.y[i]), xscale(density.x[i]));
            }
            for (let i = density.x.length - 2; i > 0; i--) {
                context.lineTo(yscale(density.y[i]), xscale(density.x[i]));
            }

            // context.closePath();

            context.fill();
            context.stroke();
            context.restore();
        }

        context.textAlign = "right";
        context.textBaseline = "middle";
        context.font = canvasFont;
        context.fillStyle = 'black';
        context.strokeStyle = 'black';

        const tickWidth = 5;
        let textWidth = size.x - tickWidth;
        context.translate(0, violinHeight * featureIndex);
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

        context.textBaseline = 'middle';

        categories.forEach((text, i) => { // features
            const pix = i * violinWidth + violinWidth / 2;
            context.save();
            // context.translate(size.x + pix, size.height);
            context.translate(size.x + pix, violinHeight + 2);
            context.rotate(-Math.PI / 2);
            context.fillText(text, 0, 0);
            context.restore();
        });
        context.textBaseline = 'top';
        context.textAlign = "middle";
        context.fillText(feature, (violinWidth * categories.length) / 2, 0);
    }

    getSize(context) {
        const {data, options} = this.props;
        const {violinHeight, violinWidth} = options;
        const categories = data.map(array => array[0].name);
        const features = data[0].map(item => item.feature);

        let maxCategoryWidth = 0;
        categories.forEach(category => {
            maxCategoryWidth = Math.max(maxCategoryWidth, context.measureText(category).width);
        });
        maxCategoryWidth += 4;

        const height = features.length * violinHeight + 4;
        const width = categories.length * violinWidth + 4;
        return {x: 30, y: maxCategoryWidth, width: width, height: height};
    }

    render() {
        this.update();
        const {saveImageEl} = this.state;
        const {data} = this.props;
        const dimension = data[0][0].dimension;
        return (<div style={{position: 'relative'}}>
            <div>
                <Typography style={{display: 'inline-block'}} component={"h4"}
                            color="textPrimary">{dimension}{this.props.subtitle &&
                <small>({this.props.subtitle})</small>}</Typography>
                <Tooltip title={"Save Image"}>
                    <IconButton aria-controls="save-image-menu" aria-haspopup="true" edge={false}
                                size={'small'}
                                aria-label="Save Image" onClick={this.handleSaveImageMenu}>
                        <PhotoCameraIcon/>
                    </IconButton>
                </Tooltip>
                <Menu
                    id="save-image-menu"
                    anchorEl={saveImageEl}
                    keepMounted
                    open={Boolean(saveImageEl)}
                    onClose={this.handleSaveImageMenuClose}
                >
                    <MenuItem onClick={e => this.handleSaveImage('png')}>PNG</MenuItem>
                    <MenuItem onClick={e => this.handleSaveImage('svg')}>SVG</MenuItem>

                </Menu>

                <Typography color="textPrimary" className="cirro-condensed" ref={this.tooltipElementRef} style={{
                    display: 'inline-block',
                    paddingLeft: 5,
                    verticalAlign: 'top',
                    whiteSpace: 'nowrap',
                    width: 500,
                    minWidth: 500,
                    maxWidth: 500,
                    textOverflow: 'ellipsis'
                }}></Typography>

            </div>
            <div ref={this.divRef}></div>
        </div>);

    }
}





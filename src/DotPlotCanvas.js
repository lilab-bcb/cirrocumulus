import {Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';
import PhotoCameraIcon from '@material-ui/icons/PhotoCamera';
import React from 'react';
import {CANVAS_FONT, SVG_FONT} from './ChartUtil';
import {drawColorScheme} from './ColorSchemeLegend';
import {numberFormat, numberFormat2f} from './formatters';
import {drawSizeLegend} from './SizeLegend';
import {stripTrailingZeros} from './util';

export const CHIP_SIZE = 12;

export function getNameWidth(array2d, context) {
    let offsets = [];
    let ncategories = 0;
    if (array2d[0].length > 0) {
        ncategories = array2d[0][0].name.length;
        for (let i = 0; i < array2d[0][0].name.length; i++) {
            offsets.push(0);
        }
    }
    array2d.forEach(array => {
        let name = array[0].name;
        for (let i = 0; i < ncategories; i++) {
            offsets[i] = Math.max(offsets[i], context.measureText(name[i]).width);
        }
    });
    for (let i = 0; i < offsets.length; i++) {
        // 4px, chip, 2px, text
        offsets[i] += 6;
        offsets[i] += CHIP_SIZE;
    }
    offsets[offsets.length - 1] += 4;
    let xTotal = 0;
    offsets.forEach(val => xTotal += val);
    return {offsets: offsets, sum: xTotal};
}


export default class DotPlotCanvas extends React.PureComponent {

    constructor(props) {
        super(props);
        this.divRef = React.createRef();
        this.tooltipElementRef = React.createRef();
        this.canvas = null;
        this.state = {saveImageEl: null};
    }


    redraw() {

        if (this.props.data == null) {
            return <div/>;
        }
        let devicePixelRatio = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            devicePixelRatio = window.devicePixelRatio;
        }

        if (this.canvas == null) {
            let onMouseMove = (event) => {
                const node = event.target;
                const maxRadius = this.props.sizeScale.range()[1];
                const rect = node.getBoundingClientRect();
                let xy = [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
                // xy[0] /= devicePixelRatio;
                // xy[1] /= devicePixelRatio;
                const col = Math.floor((xy[0] - this.size.x) / (maxRadius * 2));
                const row = Math.floor((xy[1]) / (maxRadius * 2));

                if (col >= 0 && col < this.props.data[0].length && row >= 0 && row < this.props.data.length) {
                    this.tooltipElementRef.current.innerHTML = '';
                    const array = this.props.data[row];
                    const mean = array[col].mean;
                    const fractionExpressed = array[col].fractionExpressed;
                    let meanFormatted = stripTrailingZeros(numberFormat2f(mean));
                    let percentExpressed = stripTrailingZeros(numberFormat(100 * fractionExpressed));

                    this.tooltipElementRef.current.innerHTML = 'mean: ' + meanFormatted + ', % expressed: ' + percentExpressed + ', ' + array[col].feature + ', ' + array[col].name.join(', ');
                } else {
                    this.tooltipElementRef.current.innerHTML = '';
                }
            };
            let onMouseOut = (event) => {
                this.tooltipElementRef.current.innerHTML = '';

            };
            this.canvas = document.createElement('canvas');
            this.canvas.addEventListener("mousemove", onMouseMove);
            this.canvas.addEventListener("mouseout", onMouseOut);
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
        context.font = CANVAS_FONT;

        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        this.drawContext(context, this.size);
    }

    drawContext(context, size) {
        const data2d = this.props.data;
        const colorScale = this.props.colorScale;
        const categoryColorScales = this.props.categoryColorScales;
        const sizeScale = this.props.sizeScale;
        const drawCircles = this.props.drawCircles;
        const textColor = this.props.textColor;
        const maxRadius = sizeScale.range()[1];
        let diameter = maxRadius * 2;
        // context.strokeStyle = gridColor;
        // context.lineWidth = gridThickness;

        data2d.forEach((array, j) => { // each category
            const ypix = j * diameter + (drawCircles ? maxRadius : 0);
            for (let i = 0; i < array.length; i++) { // each feature
                const mean = array[i].mean;
                const color = colorScale(mean);
                context.fillStyle = color;
                context.beginPath();
                if (drawCircles) {
                    const xpix = i * diameter + maxRadius + size.x;
                    const frac = array[i].fractionExpressed;
                    context.arc(xpix, ypix, sizeScale(frac), 0, 2 * Math.PI);
                } else {
                    const xpix = i * diameter + size.x;
                    context.rect(xpix, ypix, diameter, diameter);
                }
                context.fill();
                // context.stroke();
            }
        });

        context.textAlign = 'left';
        context.textBaseline = 'middle';

        data2d.forEach((array, i) => { // categories
            let name = array[0].name;
            const pix = i * diameter + maxRadius;
            for (let j = 0; j < name.length; j++) {
                // 4px, chip, 2px, text
                let xpixstart = size.xoffsets[j - 1] || 0;
                if (xpixstart > 0) {
                    xpixstart += 4;
                }
                const categoryColorScale = categoryColorScales[j];
                context.fillStyle = categoryColorScale(array[0].categories[j]);
                context.beginPath();
                context.rect(xpixstart, pix - maxRadius / 2 - 3, CHIP_SIZE, CHIP_SIZE);
                context.fill();
                context.stroke();
                context.fillStyle = textColor;
                context.fillText(name[j], xpixstart + 2 + CHIP_SIZE, pix);
            }
        });
        context.textAlign = 'right';
        context.textBaseline = 'top';

        data2d[0].forEach((item, i) => { // features
            const text = item.feature;
            const pix = i * diameter;
            context.save();
            context.translate(size.x + pix + 4, size.height);
            context.rotate(-Math.PI / 2);
            context.fillText(text, 0, 0);
            context.restore();
        });
        context.setTransform(1, 0, 0, 1, 0, 0);
    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    componentDidMount() {
        this.redraw();
    }

    getSize(context) {
        let maxFeatureWidth = 0;
        const array2d = this.props.data;
        array2d[0].forEach(item => {
            maxFeatureWidth = Math.max(maxFeatureWidth, context.measureText(item.feature).width);
        });
        maxFeatureWidth += 4;

        const nameWidth = getNameWidth(array2d, context);
        const maxRadius = this.props.sizeScale.range()[1];
        const diameter = maxRadius * 2;
        const height = array2d.length * diameter + 4;
        const width = array2d[0].length * diameter + 4;
        return {xoffsets: nameWidth.offsets, x: nameWidth.sum, y: maxFeatureWidth, width: width, height: height};
    }

    update() {
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        context.font = CANVAS_FONT;
        this.size = this.getSize(context);
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
            context.font = SVG_FONT;
        } else {
            canvas = document.createElement('canvas');
            context = canvas.getContext('2d');
            context.font = CANVAS_FONT;
        }

        const size = this.getSize(context);
        const colorScaleHeight = 15 + 20;
        const sizeScaleHeight = 40;
        const height = size.height + size.y + colorScaleHeight + sizeScaleHeight + 10;
        const width = Math.max(200, size.width + size.x);
        if (format === 'svg') {
            context = new window.C2S(width, height);
            context.font = SVG_FONT;
        } else {
            canvas.width = width * window.devicePixelRatio;
            canvas.height = height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.font = CANVAS_FONT;
        }
        const textColor = this.props.textColor;
        context.fillStyle = textColor === 'white' ? 'black' : 'white';
        context.fillRect(0, 0, width, height);
        this.drawContext(context, size);

        // if (format !== 'svg') {
        //     context.scale(window.devicePixelRatio, window.devicePixelRatio);
        // }

        context.translate(4, (size.height + size.y + 4));
        drawColorScheme(context, this.props.colorScale, textColor);

        context.translate(-10, (colorScaleHeight + 4));

        drawSizeLegend(context, this.props.sizeScale, 3, 150, 20, textColor);

        if (format === 'svg') {
            let svg = context.getSerializedSvg();
            // let prefix = [];
            // prefix.push('<?xml version="1.0" encoding="utf-8"?>\n');
            // prefix.push('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"' +
            //     ' "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
            // svg = prefix.join('') + svg;
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

    render() {
        this.update();
        const {saveImageEl} = this.state;
        const array2d = this.props.data;
        const dimension = array2d[0][0].dimension;

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




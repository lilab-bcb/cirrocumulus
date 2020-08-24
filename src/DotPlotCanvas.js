import {Tooltip, withStyles} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Paper from '@material-ui/core/Paper';
import Select from '@material-ui/core/Select';
import PhotoCameraIcon from '@material-ui/icons/PhotoCamera';
import {scaleLinear} from 'd3-scale';
import React from 'react';
import ColorSchemeLegend, {drawColorScheme} from './ColorSchemeLegend';
import {numberFormat} from './formatters';
import SizeLegend, {drawSizeLegend} from './SizeLegend';

const styles = theme => ({
    root: {},
    formControl: {
        display: 'block',
        margin: theme.spacing(1),
    }

});

let svgFont = '12px Helvetica,Arial,sans-serif';
let canvasFont = '12px Roboto Condensed,Helvetica,Arial,sans-serif';
let maxRadius = 9;
let minRadius = 1;
let gridColor = '#808080';
let gridThickness = 0.5;

class DotPlotCanvas extends React.PureComponent {

    constructor(props) {
        super(props);
        this.divRef = React.createRef();
        this.tooltipElementRef = React.createRef();
        this.canvas = null;
        this.state = {saveImageEl: null};
    }

    onSortOrderChanged = (event) => {
        this.props.onSortOrderChanged({name: this.props.data.name, value: event.target.value});
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
            let onMouseMove = (event) => {
                const node = event.target;
                var rect = node.getBoundingClientRect();
                let xy = [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
                xy[0] /= devicePixelRatio;
                xy[1] /= devicePixelRatio;
                const col = Math.floor((xy[0] - this.size.x) / (maxRadius * 2));
                const row = Math.floor((xy[1]) / (maxRadius * 2));

                if (col >= 0 && col < this.dotplot.values.length && row >= 0 && row < this.categories.length) {
                    this.tooltipElementRef.current.innerHTML = '';
                    let featureDatum = this.dotplot.values[col];
                    const mean = featureDatum.mean[this.categoryOrder[row]];
                    const fractionExpressed = featureDatum.fractionExpressed[this.categoryOrder[row]];
                    this.tooltipElementRef.current.innerHTML = 'mean: ' + numberFormat(mean) + ', % expressed: ' + numberFormat(100 * fractionExpressed);
                } else {
                    this.tooltipElementRef.current.innerHTML = '';
                }
            };
            let onMouseOut = (event) => {
                this.tooltipElementRef.current.innerHTML = '';

            };
            // onMouseMove = throttle(onMouseMove);
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
        context.font = canvasFont;

        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        this.drawContext(context);
    }

    drawContext(context) {
        const renamedCategories = this.props.renamedCategories || {};
        const dotplot = this.dotplot;
        const colorScale = this.colorScale;
        const features = this.features;
        const sizeScale = this.sizeScale;
        const categories = this.categories;
        const categoryOrder = this.categoryOrder;
        let diameter = maxRadius * 2;
        // context.strokeStyle = gridColor;
        // context.lineWidth = gridThickness;
        dotplot.values.forEach((datum, featureIndex) => { // each feature
            for (let i = 0; i < datum.mean.length; i++) { // each category
                const mean = datum.mean[categoryOrder[i]];
                const frac = datum.fractionExpressed[categoryOrder[i]];
                const color = colorScale(mean);
                const xpix = featureIndex * diameter + maxRadius + this.size.x;
                const ypix = i * diameter + maxRadius;
                context.fillStyle = color;
                context.beginPath();
                context.arc(xpix, ypix, sizeScale(frac), 0, 2 * Math.PI);
                context.fill();
                // context.stroke();
            }
        });
        // context.lineWidth = 1;
        context.textAlign = 'right';
        context.fillStyle = 'black';
        context.textBaseline = 'middle';

        for (let i = 0; i < categories.length; i++) {
            let category = categories[categoryOrder[i]];
            let newName = renamedCategories[category];
            if (newName != null) {
                category = newName;
            }
            const pix = i * diameter + maxRadius;
            context.fillText(category, this.size.x - 4, pix);
        }
        context.textAlign = 'right';
        context.textBaseline = 'top';
        for (let i = 0; i < features.length; i++) {
            const text = features[i];
            const pix = i * diameter;
            context.save();
            context.translate(this.size.x + pix + 4, this.size.height);
            context.rotate(-Math.PI / 2);
            context.fillText(text, 0, 0);
            context.restore();

        }

        // context.strokeStyle = gridColor;
        // context.lineWidth = gridThickness;
        //
        //
        // for (let i = 0; i < categories.length; i++) {
        //     const ypix = i * diameter;
        //     context.beginPath();
        //     context.moveTo(this.size.x + 2, ypix);
        //     context.lineTo(width, ypix);
        //     context.stroke();
        // }
        // for (let i = 0; i < features.length; i++) {
        //     const xpix = i * diameter + this.size.x + 2;
        //     context.beginPath();
        //     context.moveTo(xpix, 0);
        //     context.lineTo(xpix, dotplotHeight);
        //     context.stroke();
        // }

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
        const renamedCategories = this.props.renamedCategories || {};
        this.dotplot.values.forEach(datum => {
            maxFeatureWidth = Math.max(maxFeatureWidth, context.measureText(datum.name).width);
        });
        maxFeatureWidth += 4;
        let xoffset = 0;
        this.categories.forEach(category => {
            let renamed = renamedCategories[category];
            if (renamed != null) {
                category = renamed;
            }
            xoffset = Math.max(xoffset, context.measureText(category).width);
        });
        xoffset += 4;
        const diameter = maxRadius * 2;
        const height = this.categories.length * diameter + 4;
        const width = this.features.length * diameter + 4;
        return {x: xoffset, y: maxFeatureWidth, width: width, height: height};
    }

    update() {

        let dotplot = Object.assign({}, this.props.data);
        const renamedCategories = this.props.renamedCategories || {};
        if (dotplot != null && dotplot.selection) {
            dotplot = dotplot.selection;
        }
        this.dotplot = dotplot;
        const categories = dotplot.categories || [''];

        this.categories = categories;
        const features = dotplot.values.map(feature => feature.name);

        if (dotplot.sortBy == null) {
            dotplot.sortBy = dotplot.name;
        }
        let categoryOrder = [];
        for (let i = 0; i < categories.length; i++) {
            categoryOrder.push(i);
        }
        if (dotplot.sortBy !== dotplot.name) { // sort by feature
            let sortByDatum;
            for (let i = 0; i < dotplot.values.length; i++) {
                if (dotplot.values[i].name === dotplot.sortBy) {
                    sortByDatum = dotplot.values[i];
                    break;
                }
            }
            if (sortByDatum) {
                categoryOrder.sort((a, b) => {
                    let val1 = sortByDatum.mean[a];
                    let val2 = sortByDatum.mean[b];
                    let c = val1 === val2 ? 0 : (val1 > val2 ? -1 : 1);
                    if (c === 0) {
                        val1 = sortByDatum.fractionExpressed[a];
                        val2 = sortByDatum.fractionExpressed[b];
                        c = val1 === val2 ? 0 : (val1 > val2 ? -1 : 1);
                    }
                    return c;
                });

            }
        } else { // sort by category
            if (Object.keys(renamedCategories).length > 0) {
                categoryOrder.sort((a, b) => {
                    let val1 = categories[a];
                    let renamed1 = renamedCategories[val1];
                    if (renamed1 != null) {
                        val1 = renamed1;
                    }
                    let val2 = categories[b];
                    let renamed2 = renamedCategories[val2];
                    if (renamed2 != null) {
                        val2 = renamed2;
                    }
                    val1 = val1.toLowerCase();
                    val2 = val2.toLowerCase();
                    const num1 = parseFloat(val1);
                    const num1Valid = !isNaN(num1);
                    const num2 = parseFloat(val2);
                    const num2Valid = !isNaN(num2);

                    if (num1Valid && num2Valid) {
                        return (num1 === num2 ? 0 : (num1 < num2 ? -1 : 1));
                    }
                    if (num1Valid) {
                        return 1;
                    }
                    if (num2Valid) {
                        return -1;
                    }
                    return val1 === val2 ? 0 : (val1 > val2 ? 1 : -1);
                });
            }

        }

        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        context.font = canvasFont;
        this.features = features;
        this.size = this.getSize(context);
        this.categoryOrder = categoryOrder;
        this.colorScale = scaleLinear().domain(this.props.meanRange).range(['blue', 'red']);
        this.sizeScale = scaleLinear().domain(this.props.fractionRange).range([minRadius, maxRadius]).clamp(true);
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
        const colorScaleHeight = 40;
        const sizeScaleHeight = 40;
        const height = size.height + size.y + colorScaleHeight + sizeScaleHeight + 10;
        const width = Math.max(200, size.width + size.x);
        if (format === 'svg') {
            context = new window.C2S(width, height);
            context.font = svgFont;
        } else {
            canvas.width = width * window.devicePixelRatio;
            canvas.height = height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.fillStyle = 'white';
            context.fillRect(0, 0, width, height);
            context.font = canvasFont;
        }
        this.drawContext(context);

        if (format !== 'svg') {
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
        }
        context.translate(10, size.height + size.y + 4);
        drawColorScheme(context, 150, colorScaleHeight, this.colorScale, true);

        context.translate(-10, colorScaleHeight + 4);

        drawSizeLegend(context, this.sizeScale, 3, 150);

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
            window.saveAs(blob, this.dotplot.name + '.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, this.dotplot.name + '.png', true);
            });
        }
    };

    render() {
        this.update();
        const {saveImageEl} = this.state;
        const dotplot = this.dotplot;
        const features = this.features;
        const sortChoices = [dotplot.name].concat(features);
        return (<div style={{position: 'relative'}}>
            <Paper elevation={0}>
                <b>{dotplot.name}</b> {this.props.subtitle && <small>({this.props.subtitle})</small>}
                <Tooltip title={"Save Image"}>
                    <IconButton aria-controls="save-image-menu" aria-haspopup="true" edge={false} size={'small'}
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

                <div className="cirro-condensed" ref={this.tooltipElementRef} style={{
                    display: 'inline-block',
                    paddingLeft: 5,
                    verticalAlign: 'top',
                    whiteSpace: 'nowrap',
                    width: 500,
                    minWidth: 500,
                    maxWidth: 500,
                    textOverflow: 'ellipsis'
                }}></div>

                {this.props.onSortOrderChanged ? <FormControl className={this.props.classes.formControl}>
                    <InputLabel shrink={true}>Sort By</InputLabel>
                    <Select

                        input={<Input size={"small"}/>}
                        onChange={this.onSortOrderChanged}
                        value={dotplot.sortBy}
                    >
                        {sortChoices.map(item => (
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl> : null}
            </Paper>
            <div ref={this.divRef}></div>

            {this.props.legend && <ColorSchemeLegend style={{display: 'block', marginLeft: 10}}
                                                     width={186}
                                                     label={true} height={40}
                                                     scale={this.colorScale}/>}
            {this.props.legend && <SizeLegend style={{display: 'block'}}
                                              width={150}
                                              label={true} height={40}
                                              scale={this.sizeScale}/>}

        </div>);

    }
}

export default withStyles(styles)(DotPlotCanvas);



import {Tooltip, withStyles} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
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

class DotPlotCanvas extends React.PureComponent {

    constructor(props) {
        super(props);
        this.divRef = React.createRef();
        this.tooltipElementRef = React.createRef();
        this.backingScale = 1;
        this.maxRadius = 7;
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
        let backingScale = this.backingScale;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            backingScale = window.devicePixelRatio;
        }

        if (this.canvas == null) {
            let onMouseMove = (event) => {
                const node = event.target;
                var rect = node.getBoundingClientRect();
                let xy = [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
                xy[0] /= backingScale;
                xy[1] /= backingScale;
                const col = Math.floor((xy[0] - this.textWidth.x) / (this.maxRadius * 2));
                const row = Math.floor((xy[1]) / (this.maxRadius * 2));

                if (col >= 0 && col < this.dotplot.values.length && row >= 0 && row < this.categories.length) {
                    this.tooltipElementRef.current.innerHTML = '';
                    let featureDatum = this.dotplot.values[col];
                    const mean = featureDatum.mean[this.categoryOrder[row]];
                    const fractionExpressed = featureDatum.fractionExpressed[this.categoryOrder[row]];
                    // categories[categoryOrder[row]] + ', ' + featureDatum.name

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
        const categories = this.categories;
        const maxRadius = this.maxRadius;
        let diameter = maxRadius * 2;
        const dotplotHeight = categories.length * diameter;

        const height = dotplotHeight + this.textWidth.y;
        const width = categories.length * diameter + this.textWidth.x;
        let canvas = this.canvas;
        const context = canvas.getContext('2d');
        canvas.width = width * backingScale;
        canvas.height = height * backingScale;
        canvas.style.width = width;
        canvas.style.height = height;
        context.font = '9px Roboto Condensed,Helvetica,Arial,sans-serif';


        context
            .clearRect(0, 0, width * backingScale, height * backingScale);
        context.scale(backingScale, backingScale);
        this.drawContext(context);
    }

    drawContext(context) {
        const maxRadius = this.maxRadius;
        const dotplot = this.dotplot;
        const colorScale = this.colorScale;
        const features = this.features;
        const sizeScale = this.sizeScale;
        const categories = this.categories;
        const categoryOrder = this.categoryOrder;
        let diameter = maxRadius * 2;
        const dotplotHeight = categories.length * diameter;


        dotplot.values.forEach((datum, featureIndex) => { // each feature
            for (let i = 0; i < datum.mean.length; i++) { // each category
                const mean = datum.mean[categoryOrder[i]];
                const frac = datum.fractionExpressed[categoryOrder[i]];
                const color = colorScale(mean);
                const xpix = featureIndex * diameter + maxRadius + this.textWidth.x;
                const ypix = i * diameter + maxRadius;
                context.fillStyle = color;
                context.beginPath();
                context.arc(xpix, ypix, sizeScale(frac), 0, 2 * Math.PI);
                context.fill();
                //const text = 'mean: ' + numberFormat(mean) + ', % expressed: ' + numberFormat(100 * frac);
            }
        });

        context.textAlign = 'right';
        context.fillStyle = 'black';
        context.textBaseline = 'middle';
        for (let i = 0; i < categories.length; i++) {
            const category = categories[categoryOrder[i]];
            const pix = i * diameter + maxRadius;
            context.fillText(category, this.textWidth.x - 1, pix);
        }
        context.textAlign = 'right';
        context.textBaseline = 'top';
        for (let i = 0; i < features.length; i++) { // each category
            const text = features[i];
            const pix = i * diameter;
            context.save();
            context.translate(this.textWidth.x + pix + 2, dotplotHeight);
            context.rotate(-Math.PI / 2);
            context.fillText(text, 0, 0);
            context.restore();

        }
        // context.beginPath();
        // context.arc(xoffset, 10, 10, 0, 2 * Math.PI);
        // context.fill();
        context.setTransform(1, 0, 0, 1, 0, 0);

    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.redraw();
    }

    componentDidMount() {
        this.redraw();
    }

    getTextWidth(context) {
        let maxFeatureWidth = 0;
        this.dotplot.values.forEach(datum => {
            maxFeatureWidth = Math.max(maxFeatureWidth, context.measureText(datum.name).width);
        });
        maxFeatureWidth += 4;
        let xoffset = 0;
        this.categories.forEach(value => {
            xoffset = Math.max(xoffset, context.measureText(value).width);
        });
        xoffset += 4;
        return {x: xoffset, y: maxFeatureWidth};
    }

    update() {

        const dotplot = Object.assign({}, this.props.data);
        this.dotplot = dotplot;
        const categories = dotplot.categories || [''];
        this.categories = categories;
        const features = [];

        dotplot.values.forEach(datum => {
            features.push(datum.name);
        });
        if (dotplot.sortBy == null) {
            dotplot.sortBy = features[0];
        }
        dotplot.values = dotplot.values.filter(item => item.active);
        let colorMin = Number.MAX_VALUE;
        let colorMax = -Number.MAX_VALUE;
        let sizeMin = Number.MAX_VALUE;
        let sizeMax = -Number.MAX_VALUE;
        // set min and max values for color and size
        dotplot.values.forEach(datum => {
            datum.fractionExpressed.forEach(value => {
                sizeMin = Math.min(sizeMin, value);
                sizeMax = Math.max(sizeMax, value);
            });
            datum.mean.forEach(value => {
                colorMin = Math.min(colorMin, value);
                colorMax = Math.max(colorMax, value);
            });

        });
        let categoryOrder = [];
        for (let i = 0; i < categories.length; i++) {
            categoryOrder.push(i);
        }
        if (dotplot.sortBy !== dotplot.name) {
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
        }
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        context.font = '9px Roboto Condensed,Helvetica,Arial,sans-serif';
        this.textWidth = this.getTextWidth(context);
        if (colorMin === colorMax) {
            colorMax++;
        }
        if (sizeMin === sizeMax) {
            sizeMin = 0;
            sizeMax = 1;
        }
        this.categoryOrder = categoryOrder;

        this.features = features;

        this.colorScale = scaleLinear().domain([colorMin, colorMax]).range(['blue', 'red']);
        this.sizeScale = scaleLinear().domain([sizeMin, sizeMax]).range([1, this.maxRadius]).clamp(true);
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
        const categories = this.categories;
        const maxRadius = this.maxRadius;
        let diameter = maxRadius * 2;

        let canvas;
        if (format === 'svg') {
            context = new window.C2S(10, 10);
            context.font = '9px Helvetica,Arial,sans-serif';
        } else {
            canvas = document.createElement('canvas');
            context = canvas.getContext('2d');
            context.font = '9px Roboto Condensed,Helvetica,Arial,sans-serif';
        }
        const dotplotHeight = categories.length * diameter;
        const textWidth = this.getTextWidth(context);
        const colorScaleHeight = 40;
        const sizeScaleHeight = 40;
        const height = dotplotHeight + textWidth.y + colorScaleHeight + sizeScaleHeight + 10;
        const width = Math.max(150, categories.length * diameter + textWidth.x);
        if (format === 'svg') {
            context = new window.C2S(width, height);
            context.font = '9px Helvetica,Arial,sans-serif';
        } else {
            canvas.width = width * window.devicePixelRatio;
            canvas.height = height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.fillStyle = 'white';
            context.fillRect(0, 0, width, height);
            context.font = '9px Roboto Condensed,Helvetica,Arial,sans-serif';
        }
        this.drawContext(context);

        context.translate(10, dotplotHeight + textWidth.y + 4);

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
        const categories = this.categories;
        const sortByInputId = dotplot.name + 'sort_by';
        const sortChoices = [dotplot.name].concat(features);
        return (<div style={{position: 'relative', border: '1px solid LightGrey'}}>
            <b>{dotplot.name}</b> <small>({categories.length})</small>
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
            <FormControl className={this.props.classes.formControl}>
                <InputLabel shrink={true} htmlFor={sortByInputId}>Sort By</InputLabel>
                <Select
                    input={<Input id={sortByInputId}/>}
                    onChange={this.onSortOrderChanged}
                    value={dotplot.sortBy}
                >
                    {sortChoices.map(item => (
                        <MenuItem key={item} value={item}>{item}</MenuItem>
                    ))}
                </Select>
            </FormControl>

            <div ref={this.divRef}></div>

            <ColorSchemeLegend style={{display: 'block', marginLeft: 10}}
                               width={186}
                               label={true} height={40}
                               scale={this.colorScale}/>
            <SizeLegend style={{display: 'block'}}
                        width={150}
                        label={true} height={40}
                        scale={this.sizeScale}/>

        </div>);

    }
}

export default withStyles(styles)(DotPlotCanvas);



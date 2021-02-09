import {Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Typography from '@material-ui/core/Typography';
import PhotoCameraIcon from '@material-ui/icons/PhotoCamera';
import React from 'react';
import {CANVAS_FONT, SVG_FONT} from './ChartUtil';
import {getNameWidth} from './DotPlotCanvas';
import {numberFormat, numberFormat2f} from './formatters';
import {stripTrailingZeros} from './util';
import ViolinPlotOneFeature, {drawFeature} from './ViolinPlotOneFeature';

const yaxisWidth = 30;


export default class ViolinPlot extends React.PureComponent {

    constructor(props) {
        super(props);
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
            context.font = SVG_FONT;
        } else {
            canvas = document.createElement('canvas');
            context = canvas.getContext('2d');
            context.font = CANVAS_FONT;
        }

        const size = this.getSize(context);

        const colorScaleHeight = 15;
        const height = size.totalHeight + size.y + colorScaleHeight + 20;
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
        const textColor = 'black';
        // const textColor = this.props.textColor;
        context.fillStyle = textColor === 'white' ? 'black' : 'white';
        context.fillRect(0, 0, width, height);
        this.drawContext(context, size);

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

    drawContext(context, size) {
        const {categoryColorScales, colorScale, data, options, textColor} = this.props;
        const {violinHeight, violinWidth} = options;
        const features = data[0].map(item => item.feature);
        const categories = data.map(array => array[0].name);
        for (let i = 0; i < features.length; i++) {
            context.save();
            context.translate(0, violinHeight * i);
            drawFeature(context, size, features[i], data, colorScale, options, i === features.length - 1, categoryColorScales, textColor);
            context.textBaseline = 'top';
            context.textAlign = "middle";
            context.fillStyle = textColor;
            context.fillText(features[i], (violinWidth * categories.length) / 2, 0);
            context.restore();
        }
        context.setTransform(1, 0, 0, 1, 0, 0);
    }

    onTooltip = (item) => {
        if (item) {
            const mean = item.mean;
            const median = item.boxplotStats.median;

            const meanFormatted = stripTrailingZeros(numberFormat2f(mean));
            const medianFormatted = stripTrailingZeros(numberFormat2f(median));
            const percentExpressed = stripTrailingZeros(numberFormat(item.percentExpressed));
            this.tooltipElementRef.current.innerHTML = 'mean: ' + meanFormatted + ', median: ' + medianFormatted + ', % expressed: ' + percentExpressed + ', ' + item.name.join(', ');
        } else {
            this.tooltipElementRef.current.innerHTML = '';
        }
    };

    getSize(context) {
        const {data, options} = this.props;
        const {violinHeight, violinWidth} = options;
        const categories = data.map(array => array[0].name);
        const features = data[0].map(item => item.feature);
        const nameWidth = getNameWidth(data, context);
        const totalHeight = features.length * violinHeight + 4;
        const width = categories.length * violinWidth + 4;
        return {
            x: yaxisWidth,
            endCoordinates: nameWidth.endCoordinates,
            y: nameWidth.sum,
            width: width,
            totalHeight: totalHeight
        };
    }

    render() {

        const {saveImageEl} = this.state;
        const {categoryColorScales, colorScale, data, options, textColor} = this.props;
        const features = data[0].map(item => item.feature);
        const dimension = data[0][0].dimension;
        const dummyCanvas = document.createElement('canvas');
        const dummyContext = dummyCanvas.getContext('2d');
        dummyContext.font = CANVAS_FONT;
        const size = this.getSize(dummyContext);
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
            {features.map(feature => {
                return <ViolinPlotOneFeature onTooltip={this.onTooltip}
                                             key={feature}
                                             feature={feature} data={data}
                                             categoryColorScales={categoryColorScales}
                                             options={options}
                                             size={size}
                                             textColor={textColor}
                                             colorScale={colorScale}/>;
            })}

        </div>);

    }
}





import {Tooltip} from '@mui/material';
import IconButton from '@mui/material/IconButton';
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';
import Typography from '@mui/material/Typography';
import PhotoCameraIcon from '@mui/icons-material/PhotoCamera';
import React, {useState} from 'react';
import {CANVAS_FONT, SVG_FONT} from './ChartUtil';
import {getNameWidth} from './DotPlotCanvas';
import ViolinPlotOneFeature, {drawFeature, getViolinPlotScales} from './ViolinPlotOneFeature';

const yaxisWidth = 30;

export default function ViolinPlot(props) {

    const {categoryColorScales, colorScale, data, options, textColor} = props;
    const [saveImageEl, setSaveImageEl] = useState(null);


    function handleSaveImageMenu(event) {
        setSaveImageEl(event.currentTarget);
    }

    function handleSaveImageMenuClose(event) {
        setSaveImageEl(null);
    }

    function handleSaveImage(format) {
        setSaveImageEl(null);
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

        const size = getSize(context);

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
        // const textColor = props.textColor;
        context.fillStyle = textColor === 'white' ? 'black' : 'white';
        context.fillRect(0, 0, width, height);
        drawContext(context, size);

        if (format === 'svg') {
            let svg = context.getSerializedSvg();
            let blob = new Blob([svg], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, props.data[0][0].dimension + '.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, props.data[0][0].dimension + '.png', true);
            });
        }
    }

    function drawContext(context, size) {
        const {categoryColorScales, colorScale, data, options, textColor} = props;
        const {violinHeight, violinWidth} = options;
        const features = data[0].map(item => item.feature);
        const categories = data.map(array => array[0].name);
        for (let i = 0; i < features.length; i++) {
            context.save();
            context.translate(0, violinHeight * i);
            const {xscale, yscale} = getViolinPlotScales(data, i, options);
            drawFeature(context, size, features[i], data, colorScale, options, i === features.length - 1, categoryColorScales, textColor, xscale, yscale);
            context.textBaseline = 'top';
            context.textAlign = "middle";
            context.fillStyle = textColor;
            context.fillText(features[i], (violinWidth * categories.length) / 2, 0);
            context.restore();
        }
        context.setTransform(1, 0, 0, 1, 0, 0);
    }


    function getSize(context) {
        const {data, options} = props;
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


    const features = data[0].map(item => item.feature);
    const dimension = data[0][0].dimension;
    const dummyCanvas = document.createElement('canvas');
    const dummyContext = dummyCanvas.getContext('2d');
    dummyContext.font = CANVAS_FONT;
    const size = getSize(dummyContext);

    return (<div style={{position: 'relative'}}>
        <div>
            <Typography style={{display: 'inline-block'}} component={"h4"}
                        color="textPrimary">{dimension}{props.subtitle &&
                <small>({props.subtitle})</small>}</Typography>
            <Tooltip title={"Save Image"}>
                <IconButton aria-controls="save-image-menu" aria-haspopup="true" edge={false}
                            size={'small'}
                            aria-label="Save Image" onClick={handleSaveImageMenu}>
                    <PhotoCameraIcon/>
                </IconButton>
            </Tooltip>
            <Menu
                id="save-image-menu"
                anchorEl={saveImageEl}
                keepMounted
                open={Boolean(saveImageEl)}
                onClose={handleSaveImageMenuClose}
            >
                <MenuItem onClick={e => handleSaveImage('png')}>PNG</MenuItem>
                <MenuItem onClick={e => handleSaveImage('svg')}>SVG</MenuItem>

            </Menu>

        </div>
        {features.map(feature => {
            return <ViolinPlotOneFeature
                key={feature}
                feature={feature}
                data={data}
                categoryColorScales={categoryColorScales}
                options={options}
                size={size}
                textColor={textColor}
                colorScale={colorScale}/>;
        })}
    </div>);
}





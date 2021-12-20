import {Tooltip} from '@mui/material';
import IconButton from '@mui/material/IconButton';
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';
import Typography from '@mui/material/Typography';
import CloudDownloadIcon from '@mui/icons-material/CloudDownload';
import PhotoCameraIcon from '@mui/icons-material/PhotoCamera';
import {cumsum} from 'd3-array';
import {scaleLinear} from 'd3-scale';
import React, {useEffect, useRef, useState} from 'react';
import {CANVAS_FONT, SVG_FONT} from './ChartUtil';
import {drawColorScheme} from './ColorSchemeLegend';
import {intFormat, numberFormat, numberFormat2f} from './formatters';
import {drawSizeLegend} from './SizeLegend';
import {INTERPOLATOR_SCALING_MIN_MAX_CATEGORY, INTERPOLATOR_SCALING_MIN_MAX_FEATURE, stripTrailingZeros} from './util';
import CirroTooltip, {SCATTER_TRANSITION} from './CirroTooltip';

export const CHIP_SIZE = 12;

export function getNameWidth(array2d, context) {
    let endCoordinates = [];
    let ncategories = 0;
    if (array2d[0].length > 0) {
        ncategories = array2d[0][0].name.length;
        for (let i = 0; i < array2d[0][0].name.length; i++) {
            endCoordinates.push(0);
        }
    }
    array2d.forEach(array => {
        let name = array[0].name;
        for (let i = 0; i < ncategories; i++) {
            endCoordinates[i] = Math.max(endCoordinates[i], context.measureText(name[i]).width);
        }
    });
    for (let i = 0; i < endCoordinates.length; i++) {
        // chip, 2px, text, 4px, ...
        endCoordinates[i] += 6;
        endCoordinates[i] += CHIP_SIZE;
    }
    endCoordinates = cumsum(endCoordinates);
    endCoordinates[endCoordinates.length - 1] += 4;
    return {endCoordinates: endCoordinates, sum: endCoordinates[endCoordinates.length - 1]};
}


export default function DotPlotCanvas(props) {
    const {data, colorScale, interpolator, categoryColorScales, sizeScale, drawCircles, subtitle, textColor} = props;
    const dimension = data[0][0].dimension;
    const canvasRef = useRef();
    const [saveImageEl, setSaveImageEl] = useState(null);
    const [tip, setTip] = useState({html: ''});
    const sizeRef = useRef(null);

    useEffect(() => {
        let devicePixelRatio = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            devicePixelRatio = window.devicePixelRatio;
        }

        const canvas = canvasRef.current;
        let context = canvas.getContext('2d');

        context.font = CANVAS_FONT;
        const size = getSize(context);
        sizeRef.current = size;
        const height = size.height + size.y;
        const width = size.width + size.x;
        canvas.width = width * devicePixelRatio;
        canvas.height = height * devicePixelRatio;
        canvas.style.width = width + 'px';
        canvas.style.height = height + 'px';

        context = canvas.getContext('2d');
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        drawContext(context, size);
    });

    function exportFile() {
        const nfeatures = data.length > 0 ? data[0].length : 0;
        const ncategories = data.length;
        let text = [];
        text.push('id');
        for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
            text.push('\t');
            const name = data[categoryIndex][0].name.join('_');
            text.push(name + ':mean');
            text.push('\t');
            text.push(name + ':percent_expressed');
        }
        text.push('\n');
        for (let featureIndex = 0; featureIndex < nfeatures; featureIndex++) {
            for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
                const item = data[categoryIndex][featureIndex];
                if (categoryIndex === 0) {
                    text.push(item.feature);
                }
                text.push('\t');
                text.push(item.mean);
                text.push('\t');
                text.push(item.percentExpressed);
            }
            text.push('\n');
        }
        const blob = new Blob([text.join('')], {
            type: 'text/plain;charset=utf-8'
        });
        window.saveAs(blob, data[0][0].dimension + '.tsv');
    }

    function onMouseMove(event) {
        const node = event.target;
        const maxRadius = sizeScale.range()[1];
        const rect = node.getBoundingClientRect();
        let xy = [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
        // xy[0] /= devicePixelRatio;
        // xy[1] /= devicePixelRatio;
        const col = Math.floor((xy[0] - sizeRef.current.x) / (maxRadius * 2));
        const row = Math.floor((xy[1]) / (maxRadius * 2));

        if (col >= 0 && col < data[0].length && row >= 0 && row < data.length) {
            const array = data[row];
            const item = array[col];
            let tip = 'mean: ' + stripTrailingZeros(numberFormat2f(item.mean)) + '<br />% expressed: ' + stripTrailingZeros(numberFormat2f(item.percentExpressed)) + '<br /># cells: ' + intFormat(item.n);
            if (item.de) {
                tip += '<br />% expressed rest: ' + stripTrailingZeros(numberFormat(item.de.percentExpressed2));
                tip += '<br />log2 fold change: ' + stripTrailingZeros(numberFormat2f(item.de.foldChange));
                tip += '<br />p-value: ' + stripTrailingZeros(numberFormat2f(item.de.p));
                tip += '<br />FDR: ' + stripTrailingZeros(numberFormat2f(item.de.fdr));
                // tip += ', Mann-Whitney U : ' + stripTrailingZeros(numberFormat2f(item.results.statistic));
            }

            setTip({
                html: tip, clientX: event.clientX, clientY: event.clientY
            });
        } else {
            setTip({html: ''});
        }
    }

    function onMouseOut(event) {
        setTip({html: ''});
    }


    function drawContext(context, size) {
        const maxRadius = sizeScale.range()[1];
        const diameter = maxRadius * 2;
        // context.strokeStyle = gridColor;
        // context.lineWidth = gridThickness;
        const valueScale = interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE || interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY ? scaleLinear().range([0, 1]) : null;
        const nfeatures = data.length > 0 ? data[0].length : 0;
        const ncategories = data.length;
        let domains = null;
        if (interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
            domains = [];
            for (let featureIndex = 0; featureIndex < nfeatures; featureIndex++) {
                let min = Number.MAX_VALUE;
                let max = -Number.MAX_VALUE;
                for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
                    const array = data[categoryIndex];
                    const mean = array[featureIndex].mean;
                    min = Math.min(min, mean);
                    max = Math.max(max, mean);
                }
                domains.push([min, max]);
            }
        } else if (interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
            domains = [];
            for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
                let min = Number.MAX_VALUE;
                let max = -Number.MAX_VALUE;
                for (let featureIndex = 0; featureIndex < nfeatures; featureIndex++) {
                    const array = data[categoryIndex];
                    const mean = array[featureIndex].mean;
                    min = Math.min(min, mean);
                    max = Math.max(max, mean);
                }
                domains.push([min, max]);
            }
        }
        for (let featureIndex = 0; featureIndex < nfeatures; featureIndex++) { // each feature
            if (interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                valueScale.domain(domains[featureIndex]);
            }
            for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
                if (interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
                    valueScale.domain(domains[categoryIndex]);
                }
                const item = data[categoryIndex][featureIndex];
                let mean = item.mean;
                const ypix = categoryIndex * diameter + (drawCircles ? maxRadius : 0);
                if (valueScale) {
                    mean = valueScale(mean);
                }
                const color = colorScale(mean);
                context.fillStyle = color;
                context.beginPath();
                if (drawCircles) {
                    const xpix = featureIndex * diameter + maxRadius + size.x;
                    const frac = item.percentExpressed;
                    context.arc(xpix, ypix, sizeScale(frac), 0, 2 * Math.PI);
                } else {
                    const xpix = featureIndex * diameter + size.x;
                    context.rect(xpix, ypix, diameter, diameter);
                }
                context.fill();
            }
            // context.stroke();
        }


        context.textAlign = 'left';
        context.textBaseline = 'middle';

        data.forEach((array, i) => { // categories
            let name = array[0].name;
            const pix = i * diameter + maxRadius;
            for (let j = 0; j < name.length; j++) {
                const chipStartCoord = j === 0 ? 0 : size.endCoordinates[j - 1];
                const categoryColorScale = categoryColorScales[j];
                context.fillStyle = categoryColorScale(array[0].categories[j]);
                context.beginPath();
                context.rect(chipStartCoord, pix - maxRadius / 2 - 3, CHIP_SIZE, CHIP_SIZE);
                context.fill();
                context.stroke();
                context.fillStyle = textColor;
                context.fillText(name[j], chipStartCoord + 2 + CHIP_SIZE, pix);
            }
        });
        context.textAlign = 'right';
        context.textBaseline = 'top';

        data[0].forEach((item, i) => { // features
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


    function getSize(context) {
        let maxFeatureWidth = 0;
        data[0].forEach(item => {
            maxFeatureWidth = Math.max(maxFeatureWidth, context.measureText(item.feature).width);
        });
        maxFeatureWidth += 4;

        const nameWidth = getNameWidth(data, context);
        const maxRadius = sizeScale.range()[1];
        const diameter = maxRadius * 2;
        const height = data.length * diameter + 4;
        const width = data[0].length * diameter + 4;
        return {
            endCoordinates: nameWidth.endCoordinates, x: nameWidth.sum, y: maxFeatureWidth, width: width, height: height
        };
    }


    function handleSaveImageMenu(event) {
        setSaveImageEl(event.currentTarget);
    }

    function handleSaveImageMenuClose(event) {
        setSaveImageEl(null);
    };

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
        const colorScaleHeight = 15 + 20;
        const sizeScaleHeight = 40;
        const height = size.height + size.y + colorScaleHeight + sizeScaleHeight + 10;
        const width = Math.max(200, size.width + size.x);
        let scale = 1;
        if (format === 'svg') {
            context = new window.C2S(width, height);
            context.font = SVG_FONT;
        } else {
            canvas.width = width * window.devicePixelRatio;
            canvas.height = height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            scale = window.devicePixelRatio;
            context.scale(scale, scale);
            context.font = CANVAS_FONT;
        }

        context.fillStyle = textColor === 'white' ? 'black' : 'white';
        context.fillRect(0, 0, width, height);
        drawContext(context, size);

        // if (format !== 'svg') {
        //     context.scale(window.devicePixelRatio, window.devicePixelRatio);
        // }

        context.translate(4, scale * (size.height + size.y + 4));
        drawColorScheme(context, colorScale, textColor);
        context.translate(-10, (colorScaleHeight + 4));

        if (drawCircles) {
            drawSizeLegend(context, sizeScale, 3, 150, 20, textColor);
        }
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
            window.saveAs(blob, data[0][0].dimension + '.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, data[0][0].dimension + '.png', true);
            });
        }
    };

    return (<div style={{position: 'relative'}}>
        <div>
            <Typography style={{display: 'inline-block'}} component={"h4"}
                        color="textPrimary">{dimension}{subtitle && <small>({subtitle})</small>}</Typography>
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
            <Tooltip title={"Export"}>
                <IconButton edge={false} size={'small'} aria-label="Export" onClick={exportFile}>
                    <CloudDownloadIcon/>
                </IconButton>
            </Tooltip>
        </div>
        <div>
            <canvas ref={canvasRef} onMouseMove={onMouseMove} onMouseOut={onMouseOut}/>
            <CirroTooltip html={tip.html} clientX={tip.clientX} clientY={tip.clientY} transition={SCATTER_TRANSITION}/>
        </div>
    </div>);


}




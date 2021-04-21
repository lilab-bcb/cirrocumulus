import {Table, TableBody, TableCell, TableHead, TableRow, Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import {withStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import PhotoCameraIcon from '@material-ui/icons/PhotoCamera';
import {scaleLinear} from 'd3-scale';
import React, {useEffect, useRef, useState} from 'react';
import {CANVAS_FONT, SVG_FONT} from './ChartUtil';
import {chiSquare} from './ChiSquare';
import {fisherTest} from './FisherExact';
import {intFormat, numberFormat2f} from './formatters';

const styles = theme => ({
    table: {
        width: 'min-content',
        '& td': {padding: 6},
        '& th': {padding: 6}
    },
});

function CompositionPlot(props) {
    const [saveImageEl, setSaveImageEl] = useState(null);

    const canvasRef = useRef(null);
    const size = useRef({});
    const barWidth = 20;
    const barSpace = 10;
    const {seriesToValueToCounts, colorScale, series, textColor, title, subtitle, uniqueValues} = props;
    const barHeight = 600;

    function drawContext(context) {
        if (seriesToValueToCounts == null) {
            return;
        }
        const height = size.current.height;
        const margin = size.current.margin;
        const yScale = scaleLinear().domain([0, 1]).range([height - margin.bottom, margin.top]);
        context.textAlign = 'left';
        context.textBaseline = 'middle';
        context.textAlign = 'right';
        context.textBaseline = 'top';
        const y0 = yScale(0);
        for (let seriesIndex = 0; seriesIndex < series.length; seriesIndex++) {
            const seriesName = series[seriesIndex];
            const valueToCounts = seriesToValueToCounts[seriesName];
            let sum = 0;
            uniqueValues.forEach(uniqueValue => {
                const count = valueToCounts[uniqueValue] || 0;
                sum += count;
            });

            let yBottom = y0;
            const xPix = margin.left + barWidth * seriesIndex + (barSpace * seriesIndex);
            uniqueValues.forEach(uniqueValue => {
                const count = valueToCounts[uniqueValue];
                if (count !== undefined) {
                    let barHeight = yScale(0) - yScale(count / sum);
                    let yTop = yBottom - barHeight;
                    context.fillStyle = colorScale(uniqueValue);
                    context.fillRect(xPix, yTop, barWidth, barHeight);
                    yBottom = yTop;
                }
            });
            context.fillStyle = textColor;
            context.save();
            context.translate(xPix + 4, y0 + 4);
            context.rotate(-Math.PI / 2);
            context.fillText(seriesName, 0, 0);
            context.restore();
        }
        // yaxis
        const ticks = [0, 0.2, 0.4, 0.6, 0.8, 1];
        context.fillStyle = textColor;
        context.textBaseline = 'middle';
        context.textAlign = 'right';
        ticks.forEach(tick => {
            const pix = yScale(tick);
            context.fillText('' + 100 * tick, margin.left - 8, pix);
            context.beginPath();
            context.moveTo(margin.left - 6, pix);
            context.lineTo(margin.left - 2, pix);
            context.stroke();
        });

    }


    useEffect(() => {
        const canvas = canvasRef.current;
        let context = canvas.getContext('2d');
        context.font = CANVAS_FONT;

        let maxSeriesWidth = 0;
        series.forEach(name => maxSeriesWidth = Math.max(maxSeriesWidth, context.measureText(name).width));
        const margin = {left: 25, top: 10, bottom: maxSeriesWidth + 6, right: 4};
        size.current.margin = margin;
        const width = series.length * barWidth + series.length * barSpace + margin.left + margin.right;
        size.current.width = width;
        const height = barHeight + margin.bottom + margin.top;
        size.current.height = height;
        canvas.width = width * devicePixelRatio;
        canvas.height = height * devicePixelRatio;
        canvas.style.width = width + 'px';
        canvas.style.height = height + 'px';
        context = canvas.getContext('2d');
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, width * devicePixelRatio, height * devicePixelRatio);
        context.scale(devicePixelRatio, devicePixelRatio);
        drawContext(context);
    });

    const handleSaveImageMenu = (event) => {
        setSaveImageEl(event.currentTarget);
    };
    const handleSaveImageMenuClose = (event) => {
        setSaveImageEl(null);
    };

    const handleSaveImage = (format) => {
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
        const width = size.current.width;
        const height = size.current.height;
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
        drawContext(context, size);

        if (format === 'svg') {
            let svg = context.getSerializedSvg();
            let blob = new Blob([svg], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, 'composition.svg');
        } else {
            canvas.toBlob(blob => {
                window.saveAs(blob, 'composition.png', true);
            });
        }

    };

    const countsTable = [];

    uniqueValues.forEach(uniqueValue => {
        const counts = [];
        countsTable.push(counts);
        series.forEach(seriesName => {
            const valueToCounts = seriesToValueToCounts[seriesName];
            const count = valueToCounts[uniqueValue] || 0;
            counts.push(count);
        });
    });
    let pValue = null;
    let stat = null;
    if (countsTable.length === 2 && countsTable[0].length === 2) { // fisher exact
        pValue = fisherTest(countsTable[0][0], countsTable[0][1], countsTable[1][0], countsTable[1][1]);
        stat = 'Fisher\'s Exact';
    } else if (countsTable.length >= 2 && countsTable[0].length >= 2) {
        const result = chiSquare(countsTable);
        pValue = result.p;
        stat = 'Chi-Square';
    }

    return <>
        <div>
            <Typography style={{display: 'inline-block'}} component={"h4"}
                        color="textPrimary">{title}{subtitle &&
            <small>({subtitle})</small>}</Typography>
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
        <div style={{display: 'flex', flex: 'flex-wrap'}}>
            <canvas ref={canvasRef}></canvas>
            <div>
                <Table size={"small"} className={props.classes.table}>
                    <TableHead><TableRow><TableCell></TableCell>{series.map(item => <TableCell
                        key={item}>{item}</TableCell>)}
                    </TableRow>
                    </TableHead>
                    <TableBody>{uniqueValues.map((uniqueValue, uniqueValueIndex) => {
                        const counts = countsTable[uniqueValueIndex];
                        return <TableRow key={uniqueValue}><TableCell style={{whiteSpace: 'nowrap'}} component={"th"}>
                            <div style={{
                                display: 'inline-block',
                                width: '1em',
                                height: '1em',
                                marginRight: 2,
                                verticalAlign: 'text-bottom',
                                backgroundColor: colorScale(uniqueValue)
                            }}></div>
                            {uniqueValue}</TableCell>{series.map((seriesName, seriesIndex) => {
                            const count = counts[seriesIndex];
                            const countFormatted = intFormat(count);
                            return <TableCell key={seriesName}
                                              style={{textAlign: 'center'}}>{countFormatted}</TableCell>;
                        })}</TableRow>;
                    })}
                    </TableBody>
                </Table>
                {pValue != null &&
                <Typography color="textPrimary">{stat} p-value: {numberFormat2f(pValue)}</Typography>}
            </div>
        </div>

    </>;
}

export default withStyles(styles)(CompositionPlot);

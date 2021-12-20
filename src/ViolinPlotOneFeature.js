import Box from '@mui/material/Box';
import Typography from '@mui/material/Typography';
import {scaleLinear} from 'd3-scale';
import React, {useEffect, useRef, useState} from 'react';
import {CANVAS_FONT} from './ChartUtil';
import {CHIP_SIZE} from './DotPlotCanvas';
import {stripTrailingZeros} from './util';
import {intFormat, numberFormat2f} from './formatters';
import CirroTooltip from './CirroTooltip';


function drawAxis(context, xscale, textColor, size) {
    context.lineWidth = 1;
    const format = xscale.tickFormat(4);
    const ticks = xscale.ticks(4);
    context.textAlign = "right";
    context.textBaseline = "middle";
    context.fillStyle = textColor;
    context.strokeStyle = textColor;

    const tickWidth = 4;
    let textWidth = size.x - tickWidth;

    context.lineWidth = 0.5;
    context.beginPath();
    context.moveTo(textWidth + tickWidth, xscale(xscale.domain()[0]));
    context.lineTo(textWidth + tickWidth, xscale(xscale.domain()[1]));
    context.stroke();
    ticks.forEach(tick => {
        const pix = xscale(tick);
        context.fillText(format(tick), textWidth, pix);
        context.beginPath();
        context.moveTo(textWidth, pix);
        context.lineTo(textWidth + tickWidth, pix);
        context.stroke();
    });
    context.lineWidth = 1;
}

function Axis(props) {
    const {xscale, textColor, violinHeight, size} = props;
    const canvasRef = useRef(null);
    useEffect(() => {
        const canvas = canvasRef.current;
        const context = canvas.getContext('2d');
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, canvas.width, canvas.height);
        context.scale(window.devicePixelRatio, window.devicePixelRatio);
        drawAxis(context, xscale, textColor, size);
    });

    return (
        <canvas
            ref={canvasRef}
            width={size.x * window.devicePixelRatio}
            height={(size.y + violinHeight) * window.devicePixelRatio}
            style={{width: size.x, height: size.y + violinHeight}}
        />
    );
}

function FeatureCategory(props) {
    const {
        categoryColorScales,
        categoryIndex,
        data,
        feature,
        featureIndex,
        names,
        options,
        size,
        textColor,
        xscale,
        yscale
    } = props;
    const canvasRef = useRef(null);
    const [tip, setTip] = useState({html: ''});


    useEffect(() => {
        const canvas = canvasRef.current;
        const context = canvas.getContext('2d');
        context.setTransform(1, 0, 0, 1, 0, 0);
        context.font = CANVAS_FONT;
        context
            .clearRect(0, 0, canvas.width, canvas.height);
        context.scale(window.devicePixelRatio, window.devicePixelRatio);
        drawCategory(context, size, names, feature, featureIndex, data, categoryIndex, xscale, yscale, categoryColorScales, textColor, options, true);
    });


    return <div style={{position: 'relative', display: 'inline-block'}}>
        <canvas
            ref={canvasRef}
            onMouseOut={event => setTip({html: ''})}
            onMouseMove={event => {
                const item = data[categoryIndex][featureIndex];
                if (item) {
                    const text = 'mean: ' + stripTrailingZeros(numberFormat2f(item.mean)) + '<br /> median: ' +
                        stripTrailingZeros(numberFormat2f(item.boxplotStats.median)) + '<br /> % expressed: ' +
                        stripTrailingZeros(numberFormat2f(item.percentExpressed)) + '<br /> # cells: ' + intFormat(item.n);
                    setTip({html: text, clientX: event.clientX, clientY: event.clientY});
                } else {
                    setTip({html: ''});
                }
            }}
            width={options.violinWidth * window.devicePixelRatio}
            height={(size.y + options.violinHeight) * window.devicePixelRatio}
            style={{width: options.violinWidth, height: size.y + options.violinHeight}}
        />
        <CirroTooltip html={tip.html} clientX={tip.clientX} clientY={tip.clientY} style={{width: 150}}
                      boundsCheck={false}/>
    </div>;

}

function drawCategory(context, size, names, feature, featureIndex, data, categoryIndex, xscale, yscale, categoryColorScales, textColor, options, drawCategories) {
    const {violinScale, violinHeight, violinWidth, violinShowBoxplot} = options;
    context.save();
    context.strokeStyle = textColor;
    const lineCap = context.lineCap;
    const lineJoin = context.lineJoin;
    const item = data[categoryIndex][featureIndex];
    if (violinScale === 'width') {
        yscale = scaleLinear().domain([-item.density.max, item.density.max]).range([4, violinWidth - 4]); // horizontal position
    }
    context.lineCap = 'round';
    context.lineJoin = 'round';
    const density = item.density;
    // context.fillStyle = colorScale(item.mean);
    // context.translate(size.x + categoryIndex * violinWidth, 0);
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
    const boxplotWidth = 6;
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
    if (drawCategories) {
        // draw chips and categories
        context.textBaseline = 'middle';
        context.textAlign = "left";
        const height = violinHeight + size.y - 1;
        const name = item.name;
        const centerPix = violinWidth / 2;
        for (let j = 0; j < name.length; j++) {
            // chip, 2px, text, 4px, ...
            const chipStartCoord = j === 0 ? 0 : size.endCoordinates[j - 1];
            const categoryColorScale = categoryColorScales[j];
            const category = item.categories[j];
            context.fillStyle = categoryColorScale(category);
            context.beginPath();
            context.rect(centerPix - CHIP_SIZE + 4, height - CHIP_SIZE - chipStartCoord, CHIP_SIZE, CHIP_SIZE);
            context.fill();
            context.stroke();

            context.save();
            context.fillStyle = textColor;
            context.translate(centerPix, height - chipStartCoord - CHIP_SIZE - 2);
            context.rotate(-Math.PI / 2);
            context.fillText(name[j], 0, 0);
            context.restore();
        }
    }

}

export function getViolinPlotScales(data, featureIndex, options) {
    const {violinScale, violinHeight, violinWidth} = options;
    let xmin = Number.MAX_VALUE;
    let xmax = -Number.MAX_VALUE;
    let ymax = -Number.MAX_VALUE;
    const names = data.map(array => array[0].name); // array of arrays
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
    const xscale = scaleLinear().domain([xmin, xmax]).range([violinHeight - 1, 6]).nice(); // vertical position
    return {xscale, yscale};
}

export function drawFeature(context, size, feature, data, colorScale, options, drawCategories, categoryColorScales, textColor, xscale, yscale) {

    const features = data[0].map(item => item.feature);
    const featureIndex = features.indexOf(feature);
    const names = data.map(array => array[0].name); // array of arrays
    drawAxis(context, xscale, textColor, size);
    context.translate(size.x, 0);
    for (let categoryIndex = 0; categoryIndex < names.length; categoryIndex++) {
        context.save();
        context.translate(categoryIndex * options.violinWidth, 0);
        drawCategory(context, size, names, feature, featureIndex, data, categoryIndex, xscale, yscale, categoryColorScales, textColor, options, drawCategories);
        context.restore();
    }

    // // x axis
    // context.beginPath();
    // context.moveTo(textWidth + tickWidth, xscale(xscale.domain()[0]));
    // context.lineTo(size.width + size.x - 4, xscale(xscale.domain()[0]));
    // context.stroke();

}

export default function ViolinPlotOneFeature(props) {
    const {data, feature, textColor, size, options, onTooltip, categoryColorScales} = props;
    const features = data[0].map(item => item.feature);
    const categories = data.map(array => array[0].name);
    const featureIndex = features.indexOf(feature);
    const {xscale, yscale} = getViolinPlotScales(data, featureIndex, options);
    return (
        <Box borderColor="text.primary" border={1}
             style={{display: 'inline-block', margin: 2}}>
            <Typography color="textPrimary" component={"h4"}
                        style={{
                            marginTop: '3.2px'
                        }}>{feature}</Typography>
            <Axis xscale={xscale} textColor={textColor} violinHeight={options.violinHeight} size={size}/>
            {categories.map((category, categoryIndex) => <FeatureCategory options={options}
                                                                          key={category}
                                                                          categoryIndex={categoryIndex}
                                                                          data={data}
                                                                          feature={feature}
                                                                          featureIndex={featureIndex}
                                                                          names={categories}
                                                                          xscale={xscale}
                                                                          yscale={yscale}
                                                                          textColor={textColor}
                                                                          onTooltip={onTooltip}
                                                                          categoryColorScales={categoryColorScales}
                                                                          size={size}/>)}
        </Box>);

}






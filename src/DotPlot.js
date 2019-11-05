import {scaleLinear} from 'd3-scale';
import PropTypes from 'prop-types';
import React from 'react';
import createPlotlyComponent from 'react-plotly.js/factory';
import ColorSchemeLegend from './ColorSchemeLegend';
import {numberFormat} from './formatters';
import PlotUtil from './PlotUtil';
import SizeLegend from './SizeLegend';

const Plot = createPlotlyComponent(window.Plotly);

class DotPlot extends React.PureComponent {

    render() {
        if (this.props.data == null) {
            return <div/>;
        }
        const dotplot = this.props.data;
        let categories = dotplot.categories || [''];
        let colorMin = Number.MAX_VALUE;
        let colorMax = -Number.MAX_VALUE;
        let sizeMin = Number.MAX_VALUE;
        let sizeMax = -Number.MAX_VALUE;
        // set min and max values for color and size
        dotplot.values.forEach(datum => {
            datum.fraction_expressed.forEach(value => {
                sizeMin = Math.min(sizeMin, value);
                sizeMax = Math.max(sizeMax, value);
            });
            datum.mean.forEach(value => {
                colorMin = Math.min(colorMin, value);
                colorMax = Math.max(colorMax, value);
            });

        });


        if (colorMin === colorMax) {
            colorMax++;
        }
        if (sizeMin === sizeMax) {
            sizeMin = 0;
            sizeMax = 1;
        }
        let maxDiameter = 14;

        let colorScale = scaleLinear().domain([colorMin, colorMax]).range(['blue', 'red']);
        let sizeScale = scaleLinear().domain([sizeMin, sizeMax]).range([1, maxDiameter]).clamp(true);

        const features = [];
        dotplot.values.forEach(datum => {
            features.push(datum.name);
        });
        let traces = [];
        dotplot.values.forEach(datum => {
            const text = [];
            const color = [];
            const size = [];
            const y = [];
            for (let i = 0; i < datum.mean.length; i++) {
                y.push(datum.name);
                color.push(colorScale(datum.mean[i]));
                size.push(sizeScale(datum.fraction_expressed[i]));
                text.push('mean: ' + numberFormat(datum.mean[i]) + ', % expressed: ' + numberFormat(100 * datum.fraction_expressed[i]));
            }
            let trace = {
                x: categories,
                y: y,
                name: datum.name,
                type: 'scatter',
                text: text,
                mode: 'markers',
                sizemode: 'diameter',
                marker: {
                    color: color,
                    symbol: 'circle',
                    size: size,
                },
            };
            traces.push(trace);
        });

        let config = PlotUtil.createPlotConfig();
        let maxFeatureLength = 0;
        features.forEach(value => {
            maxFeatureLength = Math.max(maxFeatureLength, value.length);
        });

        let maxCategoryLength = 0;
        categories.forEach(value => {
            maxCategoryLength = Math.max(maxCategoryLength, value.length);
        });
        const maxFeatureWidth = 14 + maxFeatureLength * 14;
        const maxCategoryWidth = 14 + maxCategoryLength * 14;
        let layout = PlotUtil.createDotPlotLayout({
            height: 50 + maxCategoryWidth + features.length * (maxDiameter + 2),
            width: Math.max(300, 20 + maxFeatureWidth + categories.length * (maxDiameter + 2))
        });
        layout.title = {text: dotplot.name, font: {size: 12}};
        // features on y axis
        layout.margin = {l: maxFeatureWidth, b: maxCategoryWidth, t: 20, r: 0};

        return (<div style={{maxWidth: 800, overflow: 'auto', border: '1px solid LightGrey'}}>
            <Plot
                data={traces}
                layout={layout}
                config={config}
            />
            <ColorSchemeLegend style={{display: 'block'}}
                               width={300}
                               label={true} height={40}
                               scale={colorScale}/>
            <SizeLegend style={{display: 'block'}}
                        width={150}
                        label={true} height={40}
                        scale={sizeScale}/>

        </div>);
    }
}

DotPlot.propTypes = {
    data: PropTypes.object,
};


export default DotPlot;


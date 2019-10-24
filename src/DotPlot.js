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
        const data = this.props.data;
        let index = data.index || [''];
        let colorMin = Number.MAX_VALUE;
        let colorMax = -Number.MAX_VALUE;
        let sizeMin = Number.MAX_VALUE;
        let sizeMax = -Number.MAX_VALUE;
        // set min and max values
        let featureNameToValues = {};
        for (let key in data.values) {
            if (key !== 'index') {
                let values = data.values[key];
                let min;
                let max;
                let index = key.indexOf(',');
                let name = key.substring(2, index - 1);
                let type = key.substring(index + 3, key.length - 2);

                let featureValues = featureNameToValues[name];
                if (featureValues === undefined) {
                    featureValues = {};
                    featureNameToValues[name] = featureValues;
                }
                if (type === 'non_zero') {
                    min = sizeMin;
                    max = sizeMax;
                    featureValues.fraction = values;
                } else if (type === 'mean') {
                    min = colorMin;
                    max = colorMax;
                    featureValues.summary = values;
                } else {
                    console.log('Unknown type: ' + type + '.');
                }
                for (let j = 0; j < values.length; j++) {
                    min = Math.min(min, values[j]);
                    max = Math.max(max, values[j]);
                }
                if (type === 'non_zero') {
                    sizeMin = min;
                    sizeMax = max;
                } else {
                    colorMin = min;
                    colorMax = max;
                }
            }
        }
        ;
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
        let size = [];
        let color = [];
        let x = [];
        let y = [];
        let text = [];
        const features = [];

        for (let feature in featureNameToValues) {
            features.push(feature);
            let featureValues = featureNameToValues[feature];
            let summary = featureValues.summary;
            let fraction = featureValues.fraction;
            for (let j = 0; j < summary.length; j++) {
                color.push(colorScale(summary[j]));
                size.push(sizeScale(fraction[j]));
                y.push(feature);
                x.push(index[j]);
                text.push('mean: ' + numberFormat(summary[j]) + ', % non-zero: ' + numberFormat(100 * fraction[j]));
            }
        }


        let trace = {
            type: 'scatter',
            x: x,
            y: y,
            text: text,
            mode: 'markers',
            sizemode: 'diameter',
            marker: {
                color: color,
                symbol: 'circle',
                size: size,
            },
        };
        let traces = [trace];
        let config = PlotUtil.createPlotConfig();
        let maxFeatureLength = 0;
        features.forEach(value => {
            maxFeatureLength = Math.max(maxFeatureLength, value.length);
        });

        let maxCategoryLength = 0;
        index.forEach(value => {
            maxCategoryLength = Math.max(maxCategoryLength, value.length);
        });
        const maxFeatureWidth = 14 + maxFeatureLength * 14;
        const maxCategoryWidth = 14 + maxCategoryLength * 14;
        let layout = PlotUtil.createDotPlotLayout({
            height: 50 + maxCategoryWidth + features.length * (maxDiameter + 2),
            width: Math.max(300, 20 + maxFeatureWidth + index.length * (maxDiameter + 2))
        });
        layout.title = {text: data.name, font: {size: 12}};
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


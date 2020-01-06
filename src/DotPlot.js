import FormControl from '@material-ui/core/FormControl';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import withStyles from '@material-ui/core/styles/withStyles';
import {scaleLinear} from 'd3-scale';
import PropTypes from 'prop-types';
import React from 'react';
import createPlotlyComponent from 'react-plotly.js/factory';
import ColorSchemeLegend from './ColorSchemeLegend';
import {numberFormat} from './formatters';
import PlotUtil from './PlotUtil';
import SizeLegend from './SizeLegend';

const Plot = createPlotlyComponent(window.Plotly);
const styles = theme => ({
    root: {},
    formControl: {
        display: 'block',
        margin: theme.spacing(1),
    }

});

class DotPlot extends React.PureComponent {

    onSortOrderChanged = (event) => {
        this.props.onSortOrderChanged({name: this.props.data.name, value: event.target.value});
    };


    render() {
        if (this.props.data == null) {
            return <div/>;
        }
        const dotplot = Object.assign({}, this.props.data);
        dotplot.values = dotplot.values.filter(item => item.active);
        let categories = dotplot.categories || [''];
        const features = [];
        dotplot.values.forEach(datum => {
            features.push(datum.name);
        });
        if (dotplot.sortBy == null) {
            dotplot.sortBy = features[0];
        }
        let sortOrder = [];
        for (let i = 0; i < categories.length; i++) {
            sortOrder.push(i);
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
                sortOrder.sort((a, b) => {
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


        let traces = [];
        let x = [];
        for (let i = 0; i < categories.length; i++) {
            x.push(i);
        }
        dotplot.values.forEach(datum => {
            const text = [];
            const color = [];
            const size = [];
            const y = [];
            for (let i = 0; i < datum.mean.length; i++) {
                y.push(datum.name);
                const mean = datum.mean[sortOrder[i]];
                const frac = datum.fractionExpressed[sortOrder[i]];
                color.push(colorScale(mean));
                size.push(sizeScale(frac));
                text.push('mean: ' + numberFormat(mean) + ', % expressed: ' + numberFormat(100 * frac));
            }
            let trace = {
                x: x,
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

        let config = PlotUtil.createDotPlotConfig();
        let maxFeatureLength = 0;
        features.forEach(value => {
            maxFeatureLength = Math.max(maxFeatureLength, value.length);
        });

        let maxCategoryLength = 0;
        categories.forEach(value => {
            maxCategoryLength = Math.max(maxCategoryLength, value.length);
        });
        const maxFeatureWidth = maxFeatureLength * 9;
        const maxCategoryWidth = maxCategoryLength * 9;
        let layout = PlotUtil.createDotPlotLayout({
            height: 80 + maxCategoryWidth + features.length * (maxDiameter + 1),
            width: Math.max(300, maxFeatureWidth + categories.length * (maxDiameter + 1))
        });
        layout.margin = {l: maxFeatureWidth, b: maxCategoryWidth, t: 20, r: 0};

        layout.xaxis.tickvals = x;
        layout.xaxis.tickmode = 'array';
        let ticktext = [];
        for (let i = 0; i < categories.length; i++) {
            ticktext.push(categories[sortOrder[i]]);
        }
        layout.xaxis.ticktext = ticktext;
        layout.xaxis.range = [-1, categories.length];
        const sortByInputId = dotplot.name + 'sort_by';
        const sortChoices = [dotplot.name].concat(features);
        return (<div style={{maxWidth: 800, overflow: 'auto', border: '1px solid LightGrey'}}>
            <b>{dotplot.name}</b> <small>({categories.length})</small>
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

            <Plot
                data={traces}
                layout={layout}
                config={config}
            />
            <ColorSchemeLegend style={{display: 'block', marginLeft: 10}}
                               width={186}
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
    onSortOrderChanged: PropTypes.func
};


export default withStyles(styles)(DotPlot);


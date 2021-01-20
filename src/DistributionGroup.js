import {InputLabel, MenuItem, Select} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import Input from '@material-ui/core/Input';
import withStyles from '@material-ui/core/styles/withStyles';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {scaleLinear, scaleSequential} from 'd3-scale';
import {debounce} from 'lodash';
import natsort from 'natsort';
import React from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import DotPlotCanvas from './DotPlotCanvas';
import {numberFormat, numberFormat2f} from './formatters';
import {boxplotStats, density, nrd0} from './kde';
import SizeLegend from './SizeLegend';
import ViolinPlot from './ViolinPlot';

const styles = theme => ({
    root: {},
    formControl: {
        display: 'block',
        margin: theme.spacing(1),
    }
});

function updateNames(data, categoricalNames) {
    const renamedDimensions = [];
    data[0].dimensions.forEach(dimension => {
        renamedDimensions.push(categoricalNames[dimension] || {});
    });
    data.forEach(item => {
        const names = [];
        renamedDimensions.forEach((dimension, index) => {
            const nameMap = renamedDimensions[index];
            let name = item.categories[index];
            let newName = nameMap[name];
            if (newName !== undefined) {
                name = newName;
            }
            names[index] = name;
            item.name = names;
        });
    });
}

function reshapeData(data, distributionPlotOptions) {
    // create 2-d array with categories along rows and features along columns
    let categoryToFeatures = {};
    data.forEach(item => {
        let features = categoryToFeatures[item.name];
        if (features == null) {
            features = [];
            categoryToFeatures[item.name] = features;
        }
        features.push(item);
    });
    const categories = Object.keys(categoryToFeatures);
    const dimension = data[0].dimension;

    if (distributionPlotOptions.sortBy !== dimension) { // sort categories by feature
        const category2mean = {};
        const category2fractionExpressed = {};
        data.forEach(item => {
            if (item.feature === distributionPlotOptions.sortBy) {
                category2mean[item.name] = item.mean;
                category2fractionExpressed[item.name] = item.fractionExpressed;
            }
        });
        categories.sort((a, b) => {
            let val1 = category2mean[a];
            let val2 = category2mean[b];
            let c = val1 === val2 ? 0 : (val1 > val2 ? -1 : 1);
            if (c === 0) {
                val1 = category2fractionExpressed[a];
                val2 = category2fractionExpressed[b];
                c = val1 === val2 ? 0 : (val1 > val2 ? -1 : 1);
            }
            return c;
        });

    } else { // sort by category
        const sorter = natsort({insensitive: true});
        categories.sort((a, b) => {
            return sorter(a, b);
        });
    }
    let data2d = [];
    for (let i = 0; i < categories.length; i++) {
        let category = categories[i];
        data2d.push(categoryToFeatures[category]);
    }
    return data2d;
}

function getMeanAndFractionRange(result) {
    let fractionRange = [Number.MAX_VALUE, -Number.MAX_VALUE];
    let meanRange = [Number.MAX_VALUE, -Number.MAX_VALUE];
    result.forEach(feature => {
        fractionRange[0] = Math.min(feature.fractionExpressed, fractionRange[0]);
        fractionRange[1] = Math.max(feature.fractionExpressed, fractionRange[1]);
        meanRange[0] = Math.min(feature.mean, meanRange[0]);
        meanRange[1] = Math.max(feature.mean, meanRange[1]);
    });
    return {mean: meanRange, fraction: fractionRange};
}

class DistributionGroup extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {min: '', max: '', minSize: '', sizeMax: ''};
        this.updateMinSize = debounce(this.updateMinSize, 500);
        this.updateMaxSize = debounce(this.updateMaxSize, 500);
        this.updateMin = debounce(this.updateMin, 500);
        this.updateMax = debounce(this.updateMax, 500);
    }

    onMinChange = (event) => {
        this.setState({min: event.target.value});
        this.updateMin(event.target.value);
    };

    updateMin = (value) => {
        this.props.onDistributionPlotOptions({min: parseFloat(value)});
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
        this.updateMax(event.target.value);
    };

    updateMax = (value) => {
        this.props.onDistributionPlotOptions({max: parseFloat(value)});
    };


    onMinSizeChange = (event) => {
        this.setState({minSize: event.target.value});
        this.updateMinSize(event.target.value);
    };

    updateMinSize = (value) => {
        value = parseFloat(value);
        if (value > 1) {
            value /= 100; // fraction
        }
        this.props.onDistributionPlotOptions({minSize: value});
    };


    onSortOrderChanged = (event) => {
        this.props.onDistributionPlotOptions({sortBy: event.target.value});
    };

    onMaxSizeChange = (event) => {
        this.setState({sizeMax: event.target.value});
        this.updateMaxSize(event.target.value);
    };

    updateMaxSize = (value) => {
        value = parseFloat(value);
        if (value > 1) {
            value /= 100; // fraction
        }
        this.props.onDistributionPlotOptions({maxSize: value});
    };


    render() {

        const {
            textColor,
            categoryColorScales,
            distributionData,
            distributionPlotOptions,
            categoricalNames,
            selectedData,
            interpolator
        } = this.props;
        if (distributionData == null || distributionData.length === 0) {
            return null;
        }
        const meanAndFractionRange = getMeanAndFractionRange(distributionData);
        const meanRange = meanAndFractionRange.mean;
        const fractionRange = meanAndFractionRange.fraction;
        if (selectedData != null && selectedData.length > 0) {
            const selectedMeanAndFractionRange = getMeanAndFractionRange(selectedData);
            meanRange[0] = Math.min(meanRange[0], selectedMeanAndFractionRange.mean[0]);
            meanRange[1] = Math.max(meanRange[1], selectedMeanAndFractionRange.mean[1]);
            fractionRange[0] = Math.min(fractionRange[0], selectedMeanAndFractionRange.fraction[0]);
            fractionRange[1] = Math.max(fractionRange[1], selectedMeanAndFractionRange.fraction[1]);
        }
        if (meanRange[0] === meanRange[1]) {
            meanRange[1]++;
        }
        if (fractionRange[0] === fractionRange[1]) {
            fractionRange[0] = 0;
            fractionRange[1] = 1;
        }
        if (meanRange[0] > 0) {
            meanRange[0] = 0;
        }
        if (fractionRange[0] > 0) {
            fractionRange[0] = 0;
        }
        if (fractionRange[1] < 1) {
            fractionRange[1] = 1;
        }

        if (distributionPlotOptions.min != null && !isNaN(distributionPlotOptions.min)) {
            meanRange[0] = distributionPlotOptions.min;
        }
        if (distributionPlotOptions.max != null && !isNaN(distributionPlotOptions.max)) {
            meanRange[1] = distributionPlotOptions.max;
        }

        if (distributionPlotOptions.minSize != null && !isNaN(distributionPlotOptions.minSize)) {
            fractionRange[0] = distributionPlotOptions.minSize;
        }

        if (distributionPlotOptions.maxSize != null && !isNaN(distributionPlotOptions.maxSize)) {
            fractionRange[1] = distributionPlotOptions.maxSize;
        }
        if (distributionPlotOptions.sortBy == null) {
            distributionPlotOptions.sortBy = distributionData[0].dimension;
        }
        if (distributionPlotOptions.sortBy !== distributionData[0].dimension) {
            // ensure valid feature
            let found = false;
            for (let i = 0; i < distributionData.length; i++) {
                if (distributionData[i].feature === distributionPlotOptions.sortBy) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                distributionPlotOptions.sortBy = distributionData[0].dimension;
            }
        }
        const chartType = distributionPlotOptions.chartType;
        const maxRadius = 9;
        const minRadius = 1;
        const colorScale = scaleSequential(interpolator.value).domain(meanRange).clamp(true);
        const sizeScale = scaleLinear().domain(fractionRange).range([minRadius, maxRadius]).clamp(true);
        let colorMin = numberFormat(colorScale.domain()[0]);
        let colorMax = numberFormat(colorScale.domain()[1]);
        if (colorMin === colorMax) {
            colorMin = numberFormat2f(colorScale.domain()[0]);
            colorMax = numberFormat2f(colorScale.domain()[1]);
        }
        updateNames(distributionData, categoricalNames);
        if (selectedData) {
            updateNames(selectedData, categoricalNames);
        }
        const data = reshapeData(distributionData, distributionPlotOptions);
        const features = data[0].map(item => item.feature);
        if (chartType === 'violin') {
            const allData = selectedData ? distributionData.concat(selectedData) : distributionData;
            features.forEach((feature) => {
                // // ensure density is computed
                // const summary = globalFeatureSummary[feature];
                // if (summary.bandwidth == null) {
                //     summary.bandwidth = nrd0(cachedData[feature]);
                // }
                //  const bandwidth = summary.bandwidth;
                allData.forEach(item => {
                    if (item.feature === feature && item.density == null) {
                        const vector = item.vector;
                        const values = new Float32Array(vector.size());
                        for (let k = 0, n = values.length; k < n; k++) {
                            values[k] = vector.get(k);
                        }
                        item.boxplotStats = boxplotStats(values);
                        const bandwidth = nrd0(item.boxplotStats);
                        item.bandwidth = bandwidth;
                        item.density = density(values, bandwidth);
                    }
                });
            });
        }

        const sortChoices = [distributionData[0].dimension].concat(features);

        return (
            <React.Fragment>
                {chartType !== 'violin' && <DotPlotCanvas
                    categoryColorScales={categoryColorScales}
                    colorScale={colorScale}
                    sizeScale={sizeScale}
                    textColor={textColor}
                    drawCircles={chartType === 'dotplot'}
                    data={data}/>}
                {chartType === 'violin' && <ViolinPlot
                    categoryColorScales={categoryColorScales}
                    colorScale={colorScale}
                    textColor={textColor}
                    options={distributionPlotOptions}
                    textColor={textColor}
                    data={data}/>}

                {chartType !== 'violin' && selectedData != null && selectedData.length > 0 ?
                    <DotPlotCanvas
                        categoryColorScales={categoryColorScales}
                        colorScale={colorScale}
                        sizeScale={sizeScale}
                        subtitle="selection"
                        textColor={textColor}
                        drawCircles={chartType === 'dotplot'}
                        data={reshapeData(selectedData, distributionPlotOptions)}/> : null}
                {chartType === 'violin' && selectedData != null && selectedData.length > 0 ?
                    <ViolinPlot
                        categoryColorScales={categoryColorScales}
                        colorScale={colorScale}
                        subtitle="selection"
                        options={distributionPlotOptions}
                        textColor={textColor}
                        data={reshapeData(selectedData, distributionPlotOptions)}/> : null}
                {chartType !== 'violin' && <React.Fragment>
                    <ColorSchemeSelector handleInterpolator={this.props.handleInterpolator}
                                         interpolator={interpolator}/>
                    <div style={{color: textColor, width: 174}}><Typography
                        variant={"caption"}>{colorMin}</Typography><Typography
                        variant={"caption"}
                        style={{float: 'right'}}>{colorMax}</Typography></div>
                    <InputLabel shrink={true} variant={"standard"}>Custom Mean</InputLabel>
                    <TextField
                        InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                        size="small" type="text"
                        onChange={this.onMinChange} label={"Min"}
                        value={this.state.min}/>
                    <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small"
                               type="text"
                               onChange={this.onMaxChange} label={"Max"}
                               value={this.state.max}/>
                </React.Fragment>}
                {chartType === 'dotplot' && <div style={{paddingTop: 16}}>
                    <SizeLegend style={{display: 'block'}}
                                width={150}
                                textColor={textColor}
                                label={true} height={40}
                                scale={sizeScale}/>
                    <InputLabel style={{marginTop: 16}} shrink={true} variant={"standard"}>Custom Percent
                        Expressed</InputLabel>
                    <TextField InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                               size="small" type="text"
                               onChange={this.onMinSizeChange} label={"Min"}
                               value={this.state.minSize}/>
                    <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small" type="text"
                               onChange={this.onMaxSizeChange} label={"Max"}
                               value={this.state.maxSize}/>
                </div>}


                <FormControl className={this.props.classes.formControl}>
                    <InputLabel shrink={true}>Sort By</InputLabel>
                    <Select
                        input={<Input size={"small"}/>}
                        onChange={this.onSortOrderChanged}
                        value={distributionPlotOptions.sortBy}
                    >
                        {sortChoices.map(item => (
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl>

            </React.Fragment>
        );
    }

}

export default withStyles(styles)(DistributionGroup);

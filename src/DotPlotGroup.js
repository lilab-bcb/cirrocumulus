import {InputLabel} from '@material-ui/core';
import Grid from '@material-ui/core/Grid';
import Switch from '@material-ui/core/Switch';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {scaleLinear, scaleSequential} from 'd3-scale';
import {debounce} from 'lodash';
import natsort from 'natsort';
import React from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import DotPlotCanvas from './DotPlotCanvas';
import {numberFormat, numberFormat2f} from './formatters';
import SizeLegend from './SizeLegend';

function reshapeDotPlotData(data, renamedCategories, dotPlotOptions) {
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

    if (dotPlotOptions.sortBy !== dimension) { // sort categories by feature
        let category2mean = {};
        let category2fractionExpressed = {};
        data.forEach(item => {

            if (item.feature === dotPlotOptions.sortBy) {
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
        if (Object.keys(renamedCategories).length > 0) {
            const sorter = natsort();
            categories.sort((a, b) => {
                let renamed1 = renamedCategories[a];
                if (renamed1 != null) {
                    a = renamed1;
                }
                let renamed2 = renamedCategories[b];
                if (renamed2 != null) {
                    b = renamed2;
                }
                return sorter(a, b);
            });
        }
    }
    let data2d = [];
    for (let i = 0; i < categories.length; i++) {
        let category = categories[i];
        data2d.push(categoryToFeatures[category]);
    }
    return data2d;
}

function getDotPlotRange(result) {
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

export class DotPlotGroup extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {min: '', max: '', minSize: '', sizeMax: '', forceUpdate: false};
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
        this.props.onDotPlotOptions({min: parseFloat(value)});
        // this.setState({forceUpdate: !this.state.forceUpdate});
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
        this.updateMax(event.target.value);
    };

    updateMax = (value) => {
        this.props.onDotPlotOptions({max: parseFloat(value)});
        // this.setState({forceUpdate: !this.state.forceUpdate});
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
        this.props.onDotPlotOptions({minSize: value});
        // this.setState({forceUpdate: !this.state.forceUpdate});
    };


    onSortOrderChanged = (value) => {
        this.props.onDotPlotOptions({sortBy: value});
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
        this.props.onDotPlotOptions({maxSize: value});
        // this.setState({forceUpdate: !this.state.forceUpdate});
    };


    onHeatmapChange = (event) => {
        this.props.onDotPlotOptions({drawCircles: !event.target.checked});
    };


    render() {

        const {textColor, dotPlotData, dotPlotOptions, renamedCategories, selectedData, interpolator} = this.props;
        const dotPlotRange = getDotPlotRange(dotPlotData);

        const meanRange = dotPlotRange.mean;
        const fractionRange = dotPlotRange.fraction;
        if (selectedData != null && selectedData.length > 0) {
            const selectedDotPlotRange = getDotPlotRange(selectedData);
            meanRange[0] = Math.min(meanRange[0], selectedDotPlotRange.mean[0]);
            meanRange[1] = Math.max(meanRange[1], selectedDotPlotRange.mean[1]);
            fractionRange[0] = Math.min(fractionRange[0], selectedDotPlotRange.fraction[0]);
            fractionRange[1] = Math.max(fractionRange[1], selectedDotPlotRange.fraction[1]);
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

        if (dotPlotOptions.min != null && !isNaN(dotPlotOptions.min)) {
            meanRange[0] = dotPlotOptions.min;
        }
        if (dotPlotOptions.max != null && !isNaN(dotPlotOptions.max)) {
            meanRange[1] = dotPlotOptions.max;
        }

        if (dotPlotOptions.minSize != null && !isNaN(dotPlotOptions.minSize)) {
            fractionRange[0] = dotPlotOptions.minSize;
        }

        if (dotPlotOptions.maxSize != null && !isNaN(dotPlotOptions.maxSize)) {
            fractionRange[1] = dotPlotOptions.maxSize;
        }
        if (dotPlotOptions.sortBy == null) {
            dotPlotOptions.sortBy = dotPlotData[0].dimension;
        }
        if (dotPlotOptions.sortBy !== dotPlotData[0].dimension) {
            // ensure valid feature
            let found = false;
            for (let i = 0; i < dotPlotData.length; i++) {
                if (dotPlotData[i].feature === dotPlotOptions.sortBy) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                dotPlotOptions.sortBy = dotPlotData[0].dimension;
            }
        }
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
     
        return (
            <React.Fragment>
                {dotPlotData != null && dotPlotData.length > 0 ?
                    <DotPlotCanvas renamedCategories={renamedCategories}
                                   colorScale={colorScale}
                                   sizeScale={sizeScale}
                                   legend={selectedData == null || selectedData.length === 0}
                                   meanRange={meanRange}
                                   fractionRange={fractionRange}
                                   textColor={textColor}
                                   sortBy={dotPlotOptions.sortBy}
                                   drawCircles={dotPlotOptions.drawCircles}
                                   onSortOrderChanged={this.onSortOrderChanged}
                                   data={reshapeDotPlotData(dotPlotData, renamedCategories, dotPlotOptions)}/> : null}
                {selectedData != null && selectedData.length > 0 ?
                    <DotPlotCanvas renamedCategories={renamedCategories}
                                   colorScale={colorScale}
                                   sizeScale={sizeScale}
                                   subtitle="selection"
                                   legend={true}
                                   textColor={textColor}
                                   sortBy={dotPlotOptions.sortBy}
                                   drawCircles={dotPlotOptions.drawCircles}
                                   meanRange={meanRange}
                                   fractionRange={fractionRange}
                                   data={reshapeDotPlotData(selectedData, renamedCategories, dotPlotOptions)}/> : null}
                <ColorSchemeSelector handleInterpolator={this.props.handleInterpolator} interpolator={interpolator}/>
                <div style={{color: textColor, width: 174}}><Typography
                    variant={"caption"}>{colorMin}</Typography><Typography
                    variant={"caption"}
                    style={{float: 'right'}}>{colorMax}</Typography></div>
                {/*<ColorSchemeLegend style={{display: 'block', marginLeft: 10}}*/}
                {/*                   width={186}*/}
                {/*                   textColor={textColor}*/}
                {/*                   label={true} height={40}*/}
                {/*                   scale={colorScale}/>*/}

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


                {dotPlotOptions.drawCircles && <div style={{paddingTop: 16}}>
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

                <div style={{paddingTop: 16}}>
                    <Grid component="label" alignContent={"flex-start"} container alignItems="center"
                          spacing={0}>

                        <Grid item>
                            <Switch style={{display: 'inline'}} size="small" checked={!dotPlotOptions.drawCircles}
                                    onChange={this.onHeatmapChange}
                                    name="heatmap"/>
                        </Grid>
                        <Grid item><InputLabel shrink={false} variant={"standard"}>Heat Map</InputLabel></Grid>
                    </Grid>

                </div>
            </React.Fragment>
        );
    }

}

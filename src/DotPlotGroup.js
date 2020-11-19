import {InputLabel} from '@material-ui/core';
import Grid from '@material-ui/core/Grid';
import Switch from '@material-ui/core/Switch';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {scaleLinear, scaleSequential} from 'd3-scale';
import {debounce} from 'lodash';
import React from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import DotPlotCanvas from './DotPlotCanvas';
import {numberFormat} from './formatters';
import SizeLegend from './SizeLegend';

export class DotPlotGroup extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {min: '', max: '', minSize: '', sizeMax: '', forceUpdate: false, drawCircles: true};
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
        this.props.categoryItem.minCustom = parseFloat(value);
        this.setState({forceUpdate: !this.state.forceUpdate});
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
        this.updateMax(event.target.value);
    };

    updateMax = (value) => {
        this.props.categoryItem.maxCustom = parseFloat(value);
        this.setState({forceUpdate: !this.state.forceUpdate});
    };


    onMinSizeChange = (event) => {
        this.setState({minSize: event.target.value});
        this.updateMinSize(event.target.value);
    };

    updateMinSize = (value) => {
        this.props.categoryItem.minSizeCustom = parseFloat(value);
        if (this.props.categoryItem.minSizeCustom > 1) {
            this.props.categoryItem.minSizeCustom /= 100; // fraction
        }
        this.setState({forceUpdate: !this.state.forceUpdate});

    };

    onMaxSizeChange = (event) => {
        this.setState({sizeMax: event.target.value});
        this.updateMaxSize(event.target.value);
    };

    updateMaxSize = (value) => {
        this.props.categoryItem.maxSizeCustom = parseFloat(value);
        if (this.props.categoryItem.maxSizeCustom > 1) {
            this.props.categoryItem.maxSizeCustom /= 100; // fraction
        }
        this.setState({forceUpdate: !this.state.forceUpdate});
    };


    onHeatmapChange = (event) => {
        this.setState({drawCircles: !event.target.checked});
    };


    render() {

        const {textColor, categoryItem, renamedCategories, selectedData, interpolator} = this.props;
        let meanRange = categoryItem.meanRange;
        let fractionRange = categoryItem.fractionRange;

        meanRange = meanRange.slice();
        fractionRange = fractionRange.slice();
        if (selectedData != null) {

            meanRange[0] = Math.min(meanRange[0], selectedData.meanRange[0]);
            meanRange[1] = Math.max(meanRange[1], selectedData.meanRange[1]);
            fractionRange[0] = Math.min(fractionRange[0], selectedData.fractionRange[0]);
            fractionRange[1] = Math.max(fractionRange[1], selectedData.fractionRange[1]);
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


        if (categoryItem.minCustom != null && !isNaN(categoryItem.minCustom)) {
            meanRange[0] = categoryItem.minCustom;
        }
        if (categoryItem.maxCustom != null && !isNaN(categoryItem.maxCustom)) {
            meanRange[1] = categoryItem.maxCustom;
        }

        if (categoryItem.minSizeCustom != null && !isNaN(categoryItem.minSizeCustom)) {
            fractionRange[0] = categoryItem.minSizeCustom;
        }

        if (categoryItem.maxSizeCustom != null && !isNaN(categoryItem.maxSizeCustom)) {
            fractionRange[1] = categoryItem.maxSizeCustom;
        }


        const maxRadius = 9;
        const minRadius = 1;
        const colorScale = scaleSequential(interpolator.value).domain(meanRange).clamp(true);
        const sizeScale = scaleLinear().domain(fractionRange).range([minRadius, maxRadius]).clamp(true);
        return (
            <React.Fragment key={categoryItem.name}>
                <DotPlotCanvas renamedCategories={renamedCategories}
                               colorScale={colorScale}
                               sizeScale={sizeScale}
                               legend={selectedData == null}
                               meanRange={meanRange}
                               fractionRange={fractionRange}
                               textColor={textColor}
                               drawCircles={this.state.drawCircles}
                               onSortOrderChanged={this.props.onSortOrderChanged}
                               data={categoryItem}/>
                {selectedData != null ?
                    <DotPlotCanvas renamedCategories={renamedCategories}
                                   colorScale={colorScale}
                                   sizeScale={sizeScale}
                                   subtitle="selection"
                                   legend={true}
                                   textColor={textColor}
                                   drawCircles={this.state.drawCircles}
                                   meanRange={meanRange}
                                   fractionRange={fractionRange}
                                   data={selectedData}/> : null}
                <ColorSchemeSelector handleInterpolator={this.props.handleInterpolator} interpolator={interpolator}/>
                <div style={{color: textColor, width: 174}}><Typography
                    variant={"caption"}>{numberFormat(colorScale.domain()[0])}</Typography><Typography
                    variant={"caption"}
                    style={{float: 'right'}}>{numberFormat(colorScale.domain()[1])}</Typography></div>
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


                {this.state.drawCircles && <div style={{paddingTop: 16}}>
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
                            <Switch style={{display: 'inline'}} size="small" checked={!this.state.drawCircles}
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

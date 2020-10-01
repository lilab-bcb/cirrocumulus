import {InputLabel} from '@material-ui/core';
import TextField from '@material-ui/core/TextField';
import {scaleLinear, scaleSequential} from 'd3-scale';
import {interpolateReds} from 'd3-scale-chromatic';
import React from 'react';
import ColorSchemeLegend from './ColorSchemeLegend';
import DotPlotCanvas from './DotPlotCanvas';
import SizeLegend from './SizeLegend';

export class DotPlotGroup extends React.PureComponent {


    constructor(props) {
        super(props);
        this.state = {min: '', max: '', minSize: '', sizeMax: '', forceUpdate: false};
    }

    onMinChange = (event) => {
        this.setState({min: event.target.value});
    };

    onMinKeyPress = (event) => {
        if (event.key === 'Enter') {
            this.props.categoryItem.minCustom = parseFloat(event.target.value);
            this.setState({forceUpdate: !this.state.forceUpdate});
        }
    };

    onMinSizeChange = (event) => {
        this.setState({minSize: event.target.value});
    };

    onMinSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            this.props.categoryItem.minSizeCustom = parseFloat(event.target.value);
            if (this.props.categoryItem.minSizeCustom > 1) {
                this.props.categoryItem.minSizeCustom /= 100; // fraction
            }
            this.setState({forceUpdate: !this.state.forceUpdate});
        }
    };

    onMaxSizeChange = (event) => {
        this.setState({sizeMax: event.target.value});
    };

    onMaxSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            this.props.categoryItem.maxSizeCustom = parseFloat(event.target.value);
            if (this.props.categoryItem.maxSizeCustom > 1) {
                this.props.categoryItem.maxSizeCustom /= 100; // fraction
            }
            this.setState({forceUpdate: !this.state.forceUpdate});
        }
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
    };

    onMaxKeyPress = (event) => {
        if (event.key === 'Enter') {
            this.props.categoryItem.maxCustom = parseFloat(event.target.value);
            this.setState({forceUpdate: !this.state.forceUpdate});
        }
    };

    render() {

        const {textColor, categoryItem, renamedCategories, selectedData} = this.props;
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
        const colorScale = scaleSequential(interpolateReds).domain(meanRange).clamp(true);
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
                               onSortOrderChanged={this.props.onSortOrderChanged}
                               data={categoryItem}/>
                {selectedData != null ?
                    <DotPlotCanvas renamedCategories={renamedCategories}
                                   colorScale={colorScale}
                                   sizeScale={sizeScale}
                                   subtitle="selection"
                                   legend={true}
                                   textColor={textColor}
                                   meanRange={meanRange}
                                   fractionRange={fractionRange}
                                   data={selectedData}/> : null}
                <ColorSchemeLegend style={{display: 'block', marginLeft: 10}}
                                   width={186}
                                   textColor={textColor}
                                   label={true} height={40}
                                   scale={colorScale}/>
                <SizeLegend style={{display: 'block'}}
                            width={150}
                            textColor={textColor}
                            label={true} height={40}
                            scale={sizeScale}/>
                <InputLabel shrink={true} variant={"standard"}>Custom Mean</InputLabel>
                <TextField
                    InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                    size="small" type="text"
                    onKeyPress={this.onMinKeyPress}
                    onChange={this.onMinChange} label={"Min"}
                    value={this.state.min}/>
                <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small"
                           type="text"
                           onKeyPress={this.onMaxKeyPress}
                           onChange={this.onMaxChange} label={"Max"}
                           value={this.state.max}/>

                <InputLabel style={{marginTop: 16}} shrink={true} variant={"standard"}>Custom Percent
                    Expressed</InputLabel>
                <TextField InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                           size="small" type="text"
                           onKeyPress={this.onMinSizeKeyPress}
                           onChange={this.onMinSizeChange} label={"Min"}
                           value={this.state.minSize}/>
                <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small" type="text"
                           onKeyPress={this.onMaxSizeKeyPress}
                           onChange={this.onMaxSizeChange} label={"Max"}
                           value={this.state.maxSize}/>
            </React.Fragment>
        );
    }

}

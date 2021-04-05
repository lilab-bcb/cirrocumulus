import {InputLabel} from '@material-ui/core';
import Input from '@material-ui/core/Input';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';

import withStyles from '@material-ui/core/styles/withStyles';
import * as scaleChromatic from 'd3-scale-chromatic';
import React from 'react';
import ColorSchemeLegend from './ColorSchemeLegend';
import {createColorScale, fixInterpolatorName, getInterpolator, interpolators} from "./util";

const styles = theme => ({
    root: {
        display: 'flex',
        flexWrap: 'wrap',
        'flex-direction': 'column',
    },
    formControl: {
        display: 'block',
        margin: theme.spacing(1)
    },

});

function stripInterpolate(name) {
    if (name.startsWith("interpolate")) {
        name = name.substring("interpolate".length);
    }
    return name;
}


class ColorSchemeSelector extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {forceUpdate: false};
    }

    onInterpolatorChange = (event) => {
        let name = event.target.value;
        this.props.handleInterpolator(Object.assign({}, this.props.interpolator, {
            name: name,
            value: getInterpolator(name)
        }));
    };

    getScale(name) {
        return createColorScale({
            name: name,
            value: scaleChromatic[name],
            reversed: this.props.interpolator.reversed
        }).domain([0, 1]);
    }

    render() {

        const {classes, interpolator} = this.props;
        if (interpolator.reversed == null) {
            interpolator.reversed = false;
        }
        const interpolatorName = fixInterpolatorName(interpolator.name);
        const width = this.props.width || 176;
        const height = 14;
        return <>
            <InputLabel shrink={true}>Color Scheme</InputLabel>
            <Select
                input={<Input/>}
                className={classes.select}
                onChange={this.onInterpolatorChange}
                value={interpolatorName}
                multiple={false}>
                <MenuItem key="Diverging" value="Diverging" divider disabled>
                    Diverging
                </MenuItem>
                {interpolators['Diverging'].map(item => (
                    <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                        <ColorSchemeLegend width={width}
                                           label={false} height={height}
                                           scale={this.getScale(item)}/>

                    </MenuItem>))}

                <MenuItem key="Sequential (Single Hue)" value="Sequential (Single Hue)" divider disabled>
                    Sequential (Single Hue)
                </MenuItem>
                {interpolators['Sequential (Single Hue)'].map(item => (
                    <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                        <ColorSchemeLegend width={width}
                                           label={false} height={height}
                                           scale={this.getScale(item)}/>
                    </MenuItem>))}

                <MenuItem key="Sequential (Multi-Hue)" value="Sequential (Multi-Hue)" divider disabled>
                    Sequential (Multi-Hue)
                </MenuItem>
                {interpolators['Sequential (Multi-Hue)'].map(item => (
                    <MenuItem title={stripInterpolate(item)} value={item} key={item}>

                        <ColorSchemeLegend width={width}
                                           label={false} height={height}
                                           scale={this.getScale(item)}/>

                    </MenuItem>))}

                <MenuItem key="Cyclical" value="Cyclical" divider disabled>
                    Cyclical
                </MenuItem>
                {interpolators['Cyclical'].map(item => (
                    <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                        <ColorSchemeLegend width={width}
                                           label={false} height={height}
                                           scale={this.getScale(item)}/>
                    </MenuItem>))}
            </Select>
        </>;
    }
}


export default withStyles(styles)(ColorSchemeSelector);

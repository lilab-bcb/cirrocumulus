import {InputLabel} from '@mui/material';
import Input from '@mui/material/Input';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';

import withStyles from '@mui/styles/withStyles';
import * as scaleChromatic from 'd3-scale-chromatic';
import React from 'react';
import ColorSchemeLegend from './ColorSchemeLegend';
import {createColorScale, fixInterpolatorName, getInterpolator, interpolators} from "./util";

const styles = theme => ({
    root: {
        display: 'flex',
        flexWrap: 'wrap',
        'flex-direction': 'column'
    },
    formControl: {
        display: 'block',
        margin: theme.spacing(1)
    }

});

function stripInterpolate(name) {
    if (name.startsWith("interpolate")) {
        name = name.substring("interpolate".length);
    }
    return name;
}


function ColorSchemeSelector(props) {
    const {classes, interpolator, handleInterpolator, width} = props;

    function onInterpolatorChange(event) {
        let name = event.target.value;
        handleInterpolator(Object.assign({}, interpolator, {
            name: name,
            value: getInterpolator(name)
        }));
    }

    function getScale(name) {
        return createColorScale({
            name: name,
            value: scaleChromatic[name],
            reversed: interpolator == null ? false : interpolator.reversed
        }).domain([0, 1]);
    }


    if (interpolator && interpolator.reversed == null) {
        interpolator.reversed = false;
    }
    const interpolatorName = interpolator == null ? 'interpolateViridis' : fixInterpolatorName(interpolator.name);
    const _width = width || 176;
    const height = 14;
    return <>
        <InputLabel shrink={true}>Color Scheme</InputLabel>
        <Select
            disabled={interpolator == null}
            input={<Input/>}
            className={classes.select}
            onChange={onInterpolatorChange}
            value={interpolatorName}
            multiple={false}>
            <MenuItem key="Diverging" value="Diverging" divider disabled>
                Diverging
            </MenuItem>
            {interpolators['Diverging'].map(item => (
                <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                    <ColorSchemeLegend width={_width}
                                       label={false} height={height}
                                       scale={getScale(item)}/>

                </MenuItem>))}

            <MenuItem key="Sequential (Single Hue)" value="Sequential (Single Hue)" divider disabled>
                Sequential (Single Hue)
            </MenuItem>
            {interpolators['Sequential (Single Hue)'].map(item => (
                <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                    <ColorSchemeLegend width={_width}
                                       label={false} height={height}
                                       scale={getScale(item)}/>
                </MenuItem>))}

            <MenuItem key="Sequential (Multi-Hue)" value="Sequential (Multi-Hue)" divider disabled>
                Sequential (Multi-Hue)
            </MenuItem>
            {interpolators['Sequential (Multi-Hue)'].map(item => (
                <MenuItem title={stripInterpolate(item)} value={item} key={item}>

                    <ColorSchemeLegend width={_width}
                                       label={false} height={height}
                                       scale={getScale(item)}/>

                </MenuItem>))}

            <MenuItem key="Cyclical" value="Cyclical" divider disabled>
                Cyclical
            </MenuItem>
            {interpolators['Cyclical'].map(item => (
                <MenuItem title={stripInterpolate(item)} value={item} key={item}>
                    <ColorSchemeLegend width={_width}
                                       label={false} height={height}
                                       scale={getScale(item)}/>
                </MenuItem>))}
        </Select>
    </>;

}


export default withStyles(styles)(ColorSchemeSelector);

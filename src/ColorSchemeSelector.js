import Input from '@material-ui/core/Input';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';

import withStyles from '@material-ui/core/styles/withStyles';
import {scaleSequential} from 'd3-scale';
import * as scaleChromatic from 'd3-scale-chromatic';
import React from 'react';
import {connect} from 'react-redux';
import {setInterpolator} from './actions';
import ColorSchemeLegend from './ColorSchemeLegend';

const styles = theme => ({
    root: {
        display: 'flex',
        flexWrap: 'wrap',
        'flex-direction': 'column',
    },
    formControl: {
        display: 'block',
        margin: theme.spacing(1),
        minWidth: 200,
    },
    select: {
        minWidth: 200,
    },
});

let interpolators = {};
interpolators['Diverging'] = [
    'interpolateBrBG',
    'interpolatePRGn',
    'interpolatePiYG',
    'interpolatePuOr',
    'interpolateRdBu',
    'interpolateRdGy',
    'interpolateRdYlBu',
    'interpolateRdYlGn',
    'interpolateSpectral'];

interpolators['Sequential (Single Hue)'] = [
    'interpolateBlues',
    'interpolateGreens',
    'interpolateGreys',
    'interpolateOranges',
    'interpolatePurples',
    'interpolateReds'];

interpolators['Sequential (Multi-Hue)'] = [
    'interpolateViridis',
    'interpolateInferno',
    'interpolateMagma',
    'interpolatePlasma',
    'interpolateWarm',
    'interpolateCool',
    'interpolateCubehelixDefault',
    'interpolateBuGn',
    'interpolateBuPu',
    'interpolateGnBu',
    'interpolateOrRd',
    'interpolatePuBuGn',
    'interpolatePuBu',
    'interpolatePuRd',
    'interpolateRdPu',
    'interpolateYlGnBu',
    'interpolateYlGn',
    'interpolateYlOrBr',
    'interpolateYlOrRd'];

interpolators['Cyclical'] = ['interpolateRainbow', 'interpolateSinebow'];

class ColorSchemeSelector extends React.PureComponent {
    handleInterpolatorChange = (event) => {
        let name = event.target.value;
        let value = scaleChromatic[name];
        this.props.handleInterpolator({name: name, value: value});
    };

    getScale(name) {
        return scaleSequential(scaleChromatic[name]).domain([0, 1]);
    }

    render() {
        const {classes} = this.props;
        return (
            <Select
                input={<Input id="color-scheme"/>}
                className={classes.select}
                onChange={this.handleInterpolatorChange}
                value={this.props.interpolator.name}
                multiple={false}>
                <MenuItem key="Diverging" value="Diverging" divider disabled>
                    Diverging
                </MenuItem>
                {interpolators['Diverging'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={200}
                                       label={false} height={24}
                                       scale={this.getScale(item)}/>

                </MenuItem>))}

                <MenuItem key="Sequential (Single Hue)" value="Sequential (Single Hue)" divider disabled>
                    Sequential (Single Hue)
                </MenuItem>
                {interpolators['Sequential (Single Hue)'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={200}
                                       label={false} height={20}
                                       scale={this.getScale(item)}/>
                </MenuItem>))}

                <MenuItem key="Sequential (Multi-Hue)" value="Sequential (Multi-Hue)" divider disabled>
                    Sequential (Multi-Hue)
                </MenuItem>
                {interpolators['Sequential (Multi-Hue)'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={200}
                                       label={false} height={20}
                                       scale={this.getScale(item)}/>
                </MenuItem>))}

                <MenuItem key="Cyclical" value="Cyclical" divider disabled>
                    Cyclical
                </MenuItem>
                {interpolators['Cyclical'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={200}
                                       label={false} height={20}
                                       scale={this.getScale(item)}/>
                </MenuItem>))}


            </Select>
        );
    }

}

const mapStateToProps = state => {
    return {
        interpolator: state.interpolator,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleInterpolator: value => {
            dispatch(setInterpolator(value));
        },
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(ColorSchemeSelector));

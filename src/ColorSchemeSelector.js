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
import {fixInterpolatorName, getInterpolator, interpolators} from "./util";

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


class ColorSchemeSelector extends React.PureComponent {
    handleInterpolatorChange = (event) => {
        let name = event.target.value;
        this.props.handleInterpolator({name: name, value: getInterpolator(name)});
    };

    getScale(name) {
        return scaleSequential(scaleChromatic[name]).domain([0, 1]);
    }

    render() {
        const {classes} = this.props;
        let interpolator = fixInterpolatorName(this.props.interpolator.name);
        const width = 174;
        const height = 20;
        return (
            <Select
                input={<Input id="color-scheme"/>}
                className={classes.select}
                onChange={this.handleInterpolatorChange}
                value={interpolator}
                multiple={false}>
                <MenuItem key="Diverging" value="Diverging" divider disabled>
                    Diverging
                </MenuItem>
                {interpolators['Diverging'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={width}
                                       label={false} height={height}
                                       scale={this.getScale(item)}/>

                </MenuItem>))}

                <MenuItem key="Sequential (Single Hue)" value="Sequential (Single Hue)" divider disabled>
                    Sequential (Single Hue)
                </MenuItem>
                {interpolators['Sequential (Single Hue)'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={width}
                                       label={false} height={height}
                                       scale={this.getScale(item)}/>
                </MenuItem>))}

                <MenuItem key="Sequential (Multi-Hue)" value="Sequential (Multi-Hue)" divider disabled>
                    Sequential (Multi-Hue)
                </MenuItem>
                {interpolators['Sequential (Multi-Hue)'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={width}
                                       label={false} height={height}
                                       scale={this.getScale(item)}/>
                </MenuItem>))}

                <MenuItem key="Cyclical" value="Cyclical" divider disabled>
                    Cyclical
                </MenuItem>
                {interpolators['Cyclical'].map(item => (<MenuItem value={item} key={item}>
                    <ColorSchemeLegend width={width}
                                       label={false} height={height}
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

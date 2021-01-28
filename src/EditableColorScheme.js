import {InputLabel} from '@material-ui/core';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {debounce} from 'lodash';
import React from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import {numberFormat, numberFormat2f} from './formatters';


export class EditableColorScheme extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {min: '', max: ''};
        this.updateMin = debounce(this.updateMin, 500);
        this.updateMax = debounce(this.updateMax, 500);
    }

    onMinChange = (event) => {
        this.setState({min: event.target.value});
        this.updateMin(event.target.value);
    };

    updateMin = (value) => {
        this.props.onOptions({min: parseFloat(value)});
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
        this.updateMax(event.target.value);
    };

    updateMax = (value) => {
        this.props.onOptions({max: parseFloat(value)});
    };

    render() {

        const {colorScale, textColor, interpolator, onInterpolator} = this.props;
        let colorMin = numberFormat(colorScale.domain()[0]);
        let colorMax = numberFormat(colorScale.domain()[1]);
        if (colorMin === colorMax) {
            colorMin = numberFormat2f(colorScale.domain()[0]);
            colorMax = numberFormat2f(colorScale.domain()[1]);
        }
        return <React.Fragment>

            <ColorSchemeSelector handleInterpolator={onInterpolator}
                                 interpolator={interpolator}/>
            <div style={{color: textColor, width: 174}}><Typography
                variant={"caption"}>{colorMin}</Typography><Typography
                variant={"caption"}
                style={{float: 'right'}}>{colorMax}</Typography></div>
            <InputLabel shrink={true} variant={"standard"}>Custom Color Range</InputLabel>
            <TextField
                InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                size="small" type="text"
                onChange={this.onMinChange} label={"Min"}
                value={this.state.min}/>
            <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small"
                       type="text"
                       onChange={this.onMaxChange} label={"Max"}
                       value={this.state.max}/>
        </React.Fragment>;
    }
}



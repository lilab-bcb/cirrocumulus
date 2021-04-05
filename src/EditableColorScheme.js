import {InputLabel, Switch, Tooltip} from '@material-ui/core';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {debounce} from 'lodash';
import React from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import {numberFormat, numberFormat2f} from './formatters';
import {stripTrailingZeros} from './util';


export class EditableColorScheme extends React.PureComponent {

    constructor(props) {
        super(props);
        this.updateMin = debounce(this.updateMin, 500);
        this.updateMax = debounce(this.updateMax, 500);
    }

    onMinChange = (event) => {
        this.props.onMinUIChange(event.target.value);
        this.updateMin(event.target.value);
    };

    updateMin = (value) => {
        this.props.onMinChange(parseFloat(value));
    };

    onMaxChange = (event) => {
        this.props.onMaxUIChange(event.target.value);
        this.updateMax(event.target.value);
    };

    updateMax = (value) => {
        this.props.onMaxChange(parseFloat(value));
    };

    onReversedChange = (event) => {
        this.props.onInterpolator(Object.assign({}, this.props.interpolator, {reversed: event.target.checked}));
    };

    render() {
        const {domain, textColor, interpolator, onInterpolator, min, max} = this.props;
        let colorMin = "";
        let colorMax = "";
        if (domain) {
            if (!isNaN(domain[0])) {
                colorMin = stripTrailingZeros(numberFormat(domain[0]));
            }
            if (!isNaN(domain[1])) {
                colorMax = stripTrailingZeros(numberFormat(domain[1]));
            }
            if (colorMin !== '' && colorMin === colorMax) {
                colorMin = stripTrailingZeros(numberFormat2f(domain[0]));
                colorMax = stripTrailingZeros(numberFormat2f(domain[1]));
            }
        }
        const width = 176;
        return <>
            <ColorSchemeSelector handleInterpolator={onInterpolator}
                                 interpolator={interpolator}/>
            <>
                <div style={{color: textColor, width: width}}><Typography
                    variant={"caption"}>{colorMin}</Typography><Typography
                    variant={"caption"}
                    style={{float: 'right'}}>{colorMax}</Typography></div>
                <InputLabel disabled={domain == null} shrink={true} variant={"standard"}>Custom Color
                    Range</InputLabel>
                <TextField
                    InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                    size="small" type="text"
                    disabled={domain == null}
                    onChange={this.onMinChange} label={"Min"}
                    value={min}/>
                <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small"
                           type="text"
                           disabled={domain == null}
                           onChange={this.onMaxChange} label={"Max"}
                           value={max}/>
            </>
            <Tooltip title={"Select to invert the color order"}>
                <div><FormControlLabel
                    control={
                        <Switch
                            checked={interpolator.reversed}
                            onChange={this.onReversedChange}
                        />
                    }
                    label="Reverse Colors"
                /></div>
            </Tooltip>
        </>;
    }
}



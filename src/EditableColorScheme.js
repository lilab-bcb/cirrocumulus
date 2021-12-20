import {InputLabel, Switch, Tooltip} from '@mui/material';
import FormControlLabel from '@mui/material/FormControlLabel';
import TextField from '@mui/material/TextField';
import Typography from '@mui/material/Typography';
import {debounce} from 'lodash';
import React, {useEffect, useMemo} from 'react';
import ColorSchemeSelector from './ColorSchemeSelector';
import {numberFormat, numberFormat2f} from './formatters';
import {stripTrailingZeros} from './util';


export function EditableColorScheme(props) {
    const {
        domain,
        textColor,
        interpolator,
        onMinChange,
        onMaxChange,
        onMinUIChange,
        onMaxUIChange,
        onInterpolator,
        min,
        max
    } = props;


    function updateMin(value) {
        onMinChange(parseFloat(value));
    }

    function updateMax(value) {
        onMaxChange(parseFloat(value));
    }

    function onReversedChange(event) {
        onInterpolator(Object.assign({}, interpolator, {reversed: event.target.checked}));
    }

    const updateMinDebounced = useMemo(() => debounce(updateMin, 5000), []);
    const updateMaxDebounced = useMemo(() => debounce(updateMax, 500), []);

    useEffect(() => {
        return () => {
            updateMinDebounced.cancel();
            updateMaxDebounced.cancel();
        };
    }, [updateMinDebounced, updateMaxDebounced]);

    function handleMin(event) {
        onMinUIChange(event.target.value);
        updateMinDebounced(event.target.value);
    }

    function handleMax(event) {
        onMaxUIChange(event.target.value);
        updateMaxDebounced(event.target.value);
    }


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
                onChange={handleMin} label={"Min"}
                value={min}/>
            <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small"
                       type="text"
                       disabled={domain == null}
                       onChange={handleMax} label={"Max"}
                       value={max}/>
        </>
        <Tooltip title={"Select to invert the color order"}>
            <div><FormControlLabel
                control={
                    <Switch
                        disabled={domain == null}
                        checked={interpolator == null ? false : interpolator.reversed}
                        onChange={onReversedChange}
                    />
                }
                label="Reverse Colors"
            /></div>
        </Tooltip>
    </>;

}



import {InputLabel, Switch} from '@mui/material';
import FormControlLabel from '@mui/material/FormControlLabel';
import TextField from '@mui/material/TextField';
import {debounce} from 'lodash';
import React, {useEffect, useMemo, useState} from 'react';
import SizeLegend from './SizeLegend';


export function EditableSizeLegend(props) {
    const {onOptions, sizeScale, reversed, onReversedChange, showReversed, textColor} = props;

    const [minSize, setMinSize] = useState('');
    const [maxSize, setMaxSize] = useState('');


    function updateMinSize(value) {
        value = parseFloat(value);
        onOptions({minSize: value});
    }


    function updateMaxSize(value) {
        value = parseFloat(value);
        onOptions({maxSize: value});
    }

    function handleReversedChange(event) {
        onReversedChange(event.target.checked);
    }

    const updateMinSizeDebounced = useMemo(() => debounce(updateMinSize, 500), []);
    const updateMaxSizeDebounced = useMemo(() => debounce(updateMaxSize, 500), []);
    useEffect(() => {
        return () => {
            updateMinSizeDebounced.cancel();
            updateMaxSizeDebounced.cancel();
        };
    }, [updateMinSizeDebounced, updateMaxSizeDebounced]);

    function onMinSizeChange(event) {
        setMinSize(event.target.value);
        updateMinSizeDebounced(event.target.value);
    }

    function onMaxSizeChange(event) {
        setMaxSize(event.target.value);
        updateMaxSizeDebounced(event.target.value);
    }


    return <>
        <SizeLegend style={{display: 'block'}}
                    width={174}
                    textColor={textColor}
                    label={true} height={40}
                    scale={sizeScale}/>
        {showReversed && <div><FormControlLabel
            control={
                <Switch
                    checked={reversed}
                    onChange={handleReversedChange}
                />
            }
            label="Reverse Sizes"
        /></div>}
        <InputLabel style={{marginTop: 16}} shrink={true} variant={"standard"}>Custom Size Range</InputLabel>
        <TextField InputLabelProps={{shrink: true}} style={{width: 90, marginRight: 4}}
                   size="small" type="text"
                   onChange={onMinSizeChange} label={"Min"}
                   value={minSize}/>
        <TextField InputLabelProps={{shrink: true}} style={{width: 90}} size="small" type="text"
                   onChange={onMaxSizeChange} label={"Max"}
                   value={maxSize}/>
    </>;
}



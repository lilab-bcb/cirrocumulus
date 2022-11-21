import {Switch, Typography} from '@mui/material';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';
import TextField from '@mui/material/TextField';
import {debounce} from 'lodash';
import React, {useMemo} from 'react';
import FormControlLabel from '@mui/material/FormControlLabel';

export default function MeasureFilter(props) {
  const {handleUpdate, datasetFilter, name} = props;
  let filters = datasetFilter[name];
  if (filters == null) {
    // default
    filters = {
      invert: false,
      operation: ['>', '<'],
      value: [NaN, NaN],
      uiValue: ['', ''],
    };
    datasetFilter[name] = filters;
  }

  function handleValueUpdate(index) {
    const priorValue = filters.value[index];
    filters.value[index] = parseFloat(filters.uiValue[index]);

    handleUpdate({
      name: name,
      filter: filters,
      update: priorValue !== filters.value[index],
    });
  }

  const handleValueUpdateDebouncedFunc = useMemo(
    () => debounce(handleValueUpdate, 500),
    [name],
  );

  function handleOperationChanged(event, index) {
    filters.operation[index] = event.target.value;
    handleUpdate({
      name: name,
      filter: filters,
      update: true,
    });
  }

  function handleValueChange(event, index) {
    filters.uiValue[index] = event.target.value;
    handleUpdate({
      name: name,
      filter: filters,
      update: false,
    });
    handleValueUpdateDebouncedFunc(index);
  }

  function handleInvertChange(event) {
    filters.invert = event.target.checked;
    handleUpdate({
      name: name,
      filter: filters,
      update: true,
    });
  }

  return (
    <div>
      <Typography
        gutterBottom={false}
        component={'h2'}
        style={{textTransform: 'uppercase', letterSpacing: '0.1em'}}
      >
        Filter
      </Typography>
      {[0, 1].map((index) => {
        return (
          <div key={index}>
            <Select
              autoWidth
              size={'small'}
              style={{marginRight: 6}}
              value={filters.operation[index]}
              onChange={(event) => handleOperationChanged(event, index)}
            >
              <MenuItem value={''}></MenuItem>
              <MenuItem value={'>'}>{'>'}</MenuItem>
              <MenuItem value={'<'}>{'<'}</MenuItem>
              <MenuItem value={'='}>{'='}</MenuItem>
              <MenuItem value={'>='}>{'>='}</MenuItem>
              <MenuItem value={'<='}>{'<='}</MenuItem>
              <MenuItem value={'!='}>{'!='}</MenuItem>
            </Select>

            <TextField
              autoComplete={'off'}
              size={'small'}
              onChange={(event) => handleValueChange(event, index)}
              value={filters.uiValue[index]}
              style={{width: 60}}
            />
          </div>
        );
      })}
      <div style={{alignItems: 'flex-end'}}>
        <FormControlLabel
          control={<Switch size={'small'} onChange={handleInvertChange} />}
          label="Invert Filter"
        />
      </div>
    </div>
  );
}

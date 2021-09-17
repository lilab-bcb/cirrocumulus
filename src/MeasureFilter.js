import {InputLabel} from '@mui/material';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';
import TextField from '@mui/material/TextField';
import {debounce} from 'lodash';
import React from 'react';


class MeasureFilter extends React.PureComponent {

    constructor(props) {
        super(props);
        this.handleValueUpdate = debounce(this.handleValueUpdate, 500);
    }

    getFilter() {
        let filter = this.props.datasetFilter[this.props.name];
        if (filter == null) {
            filter = {operation: ">", value: NaN, uiValue: ""};
            this.props.datasetFilter[this.props.name] = filter;
        }
        return filter;
    }

    handleValueUpdate = () => {
        const filter = this.getFilter();
        let value = parseFloat(filter.uiValue);
        this.props.handleUpdate({name: this.props.name, value: value, update: true});

    };

    handleOperationChanged = (event) => {
        const operation = event.target.value;
        this.props.handleUpdate({name: this.props.name, operation: operation, update: true});
    };

    handleValueChange = (event) => {
        const filter = this.getFilter();
        filter.uiValue = event.target.value;
        this.props.handleUpdate({
            name: this.props.name,
            operation: filter.operation,
            value: filter.value,
            update: false
        });
        this.handleValueUpdate();
    };

    render() {
        const {name} = this.props;
        const filter = this.getFilter();
        const id = name + '_filter';

        return (

            <div style={{display: 'flex'}}>
                <InputLabel id={id + '_label'}>Filter</InputLabel>
                <Select
                    label={"Filter"}
                    size={"small"}
                    labelId={id + '_label'}
                    id={id}
                    style={{marginRight: 6}}
                    value={filter.operation}
                    onChange={this.handleOperationChanged}
                >
                    <MenuItem value={""}></MenuItem>
                    <MenuItem value={">"}>{">"}</MenuItem>
                    <MenuItem value={"<"}>{"<"}</MenuItem>
                    <MenuItem value={"="}>{"="}</MenuItem>
                    <MenuItem value={">="}>{">="}</MenuItem>
                    <MenuItem value={"<="}>{"<="}</MenuItem>
                    <MenuItem value={"!="}>{"!="}</MenuItem>
                </Select>

                <TextField size={"small"}
                           onChange={this.handleValueChange} value={filter.uiValue} style={{maxWidth: 60}}/>

            </div>
        );
    }
}

export default MeasureFilter;




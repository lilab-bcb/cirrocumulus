import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import TextField from '@material-ui/core/TextField';
import React from 'react';


class MeasureFilter extends React.PureComponent {

    getFilter() {
        let filter = this.props.datasetFilter[this.props.name];
        if (filter == null) {
            filter = {operation: ">", value: NaN, uiValue: ""};
            this.props.datasetFilter[this.props.name] = filter;
        }
        return filter;
    }

    handleOperationChanged = (event) => {
        const operation = event.target.value;
        this.props.handleUpdate({name: this.props.name, operation: operation, update: true});
    };


    handleValueKeyPress = (event) => {
        if (event.key === 'Enter') {
            const filter = this.getFilter();
            let value = parseFloat(filter.uiValue);
            this.props.handleUpdate({name: this.props.name, value: value, update: true});
        }
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
    };

    render() {
        const {name} = this.props;
        const filter = this.getFilter();
        const id = name + '_filter';

        return (

            <div style={{display: 'inline-flex', paddingLeft: 10}}>

                <Select
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

                <TextField onKeyPress={this.handleValueKeyPress}
                           onChange={this.handleValueChange} value={filter.uiValue} style={{maxWidth: 90}}/>
            </div>
        );
    }
}

export default MeasureFilter;




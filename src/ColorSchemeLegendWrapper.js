import {InputLabel} from '@material-ui/core';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import ContinuousLegend from './ContinuousLegend';
import MeasureFilter from './MeasureFilter';

class ColorSchemeLegendWrapper extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {min: '', max: ''};
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.name !== this.props.name) {
            this.setState({
                min: '',
                max: ''
            });
        }
    }

    onMinChange = (event) => {
        this.setState({min: event.target.value});
    };

    onMinKeyPress = (event) => {
        if (event.key === 'Enter') {
            let value = parseFloat(event.target.value);
            if (isNaN(value)) {
                value = this.props.globalFeatureSummary[this.props.name].min;
            }
            let domain = this.props.scale.domain();
            domain[0] = value;
            this.setState({min: event.target.value});
            this.props.handleDomain({name: this.props.name, domain: domain});

        }
    };

    onMaxChange = (event) => {
        this.setState({max: event.target.value});
    };

    onMaxKeyPress = (event) => {
        if (event.key === 'Enter') {
            let value = parseFloat(event.target.value);
            if (isNaN(value)) {
                value = this.props.globalFeatureSummary[this.props.name].max;
            }
            let domain = this.props.scale.domain();
            domain[1] = value;
            this.setState({max: event.target.value});
            this.props.handleDomain({name: this.props.name, domain: domain});

        }
    };

    render() {
        const {scale, name, nObs, nObsSelected, featureSummary, globalFeatureSummary, maxHeight, datasetFilter, handleUpdate, selected} = this.props;
        let style = {display: 'inline-block', verticalAlign: 'top'};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }

        return (
            <div className="cirro-condensed" style={style}>
                {/*ContinuousLegend shows stats table */}
                <ContinuousLegend name={name}
                                  selected={selected}
                                  summary={scale.summary}
                                  featureSummary={featureSummary}
                                  nObs={nObs}
                                  nObsSelected={nObsSelected}
                                  globalFeatureSummary={globalFeatureSummary}
                                  maxHeight={maxHeight}></ContinuousLegend>
                {name !== '__count' && this.props.handleDomain &&
                <InputLabel shrink={true} variant={"standard"}>Custom Color Scale</InputLabel>}
                {name !== '__count' && this.props.handleDomain &&
                <TextField margin="none" style={{width: 90, marginRight: 4}} size="small" type="text"
                           onKeyPress={this.onMinKeyPress}
                           onChange={this.onMinChange} placeholder="Min"
                           value={this.state.min}/>}
                {name !== '__count' && this.props.handleDomain &&
                <TextField margin="none" style={{width: 90}} size="small" type="text" onKeyPress={this.onMaxKeyPress}
                           onChange={this.onMaxChange} placeholder="Max"
                           value={this.state.max}/>}
                {name !== '__count' &&
                <MeasureFilter datasetFilter={datasetFilter} name={name} handleUpdate={handleUpdate}/>}
                {/*{this.props.showColorScheme &&*/}
                {/*<ColorSchemeLegend width={this.props.width} height={this.props.height} style={{paddingLeft: 10}}*/}
                {/*                   scale={scale}*/}
                {/*                   label={this.props.label}></ColorSchemeLegend>}*/}
            </div>);

    }
}

export default ColorSchemeLegendWrapper;

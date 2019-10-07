import {Switch} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import ListItemText from '@material-ui/core/ListItemText';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';

import withStyles from '@material-ui/core/styles/withStyles';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {
    setBinSummary,
    setBinValues,
    setFeatures,
    setMarkerOpacity,
    setMarkerOpacityUI,
    setMarkerSize,
    setMarkerSizeUI,
    setNumberOfBins,
    setNumberOfBinsUI,
    setView3d,
    setViewName,
} from './actions';

import Autocomplete from './Autocomplete';
import ColorSchemeSelector from './ColorSchemeSelector';

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

class EmbedForm extends React.PureComponent {

    onMarkerSizeChange = (event) => {
        this.props.handleMarkerSizeUI(event.target.value);
    };

    onMarkerSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            let markerSize = parseFloat(event.target.value);
            if (markerSize > 0) {
                this.props.handleMarkerSize(markerSize);
            }
        }
    };

    onMarkerOpacityChange = (event) => {
        this.props.handleMarkerOpacityUI(event.target.value);
    };

    onMarkerOpacityKeyPress = (event) => {
        if (event.key === 'Enter') {
            let opacity = parseFloat(event.target.value);
            if (opacity > 0 && opacity <= 1) {
                this.props.handleMarkerOpacity(opacity);
            }
        }
    };

    onNumberOfBinsChange = (event) => {
        this.props.handleNumberOfBinsUI(event.target.value);
    };

    onNumberOfBinsKeyPress = (event) => {
        if (event.key === 'Enter') {
            let numberOfBins = parseInt(event.target.value);
            if (numberOfBins > 0) {
                this.props.handleNumberOfBins(numberOfBins);
            }
        }
    };

    onBinSummaryChange = (event) => {
        this.props.handleBinSummary(event.target.value);
    };

    handleBinValuesChange = (event) => {
        this.props.handleBinValues(event.target.checked);
    };

    handleViewChange = (event) => {
        const viewName = event.target.value;
        this.props.handleViewName(viewName);
    };

    handle3dChange = (event) => {
        const is3d = event.target.checked;
        this.props.handleView3d(is3d);
    };


    render() {
        const {classes, selectedFeatures, selectedGroupBy, viewName, view3d, numberOfBins, markerSize, markerOpacity, binValues, binSummary, features, views, obsCat, obs} = this.props;
        let selectedView = {dimensions: 2};
        for (let i = 0; i < views.length; i++) {
            if (views[i].name === viewName) {
                selectedView = views[i];
                break;
            }
        }


        let is3dAvailable = selectedView.dimensions > 2;
        const summaryOptions = [
            {value: 'min', label: 'Minimum'},
            {value: 'max', label: 'Maximum'},
            {value: 'mean', label: 'Mean'},
            {value: 'sum', label: 'Sum'}];
        let featureValue = selectedFeatures.concat(selectedGroupBy);
        featureValue = featureValue.map(item => {
            return {label: item, value: item};
        });
        let metadataOptions = obs.map(item => {
            return {label: item, value: item};
        });
        metadataOptions = metadataOptions.concat(obsCat.map(item => {
            return {label: item, value: item, categorical: true};
        }));
        metadataOptions.sort((a, b) => {
            a = a.label.toLowerCase();
            b = b.label.toLowerCase();
            return a < b ? -1 : (a == b ? 0 : 1);
        });

        let allOptions = [{label: 'Annotations', options: metadataOptions}, {
            label: 'Variables',
            options: features.map(item => {
                return {label: item, value: item};
            })
        }];
        let defaultOptions = [{
            label: 'Annotations', options: metadataOptions
        }, {
            label: 'Variables',
            options: [{isDisabled: true, label: 'Type to search', value: ''}]
        }]
        return (
            <div className={classes.root}>
                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="view">Embedding</InputLabel>
                    <Select
                        className={classes.select}
                        input={<Input id="view"/>}
                        onChange={this.handleViewChange}
                        value={viewName}
                        multiple={false}>
                        {views.map(view => (
                            <MenuItem key={view.name} value={view.name}>
                                <ListItemText primary={view.name}/>
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>


                <FormControl className={classes.formControl}>
                    <Autocomplete label="Features" options={allOptions}
                                  defaultOptions={defaultOptions} value={featureValue}
                                  onChange={this.props.handleFeatures}
                                  isMulti={true}/>
                </FormControl>


                {is3dAvailable && <FormControlLabel
                    disabled={!is3dAvailable}
                    control={
                        <Switch
                            checked={view3d}
                            onChange={this.handle3dChange}
                        />
                    }
                    label="3-d"
                />}


                <TextField type="number" step="2" min="0.1" max="30" onKeyPress={this.onMarkerSizeKeyPress}
                           onChange={this.onMarkerSizeChange} label="Marker Size"
                           className={classes.formControl} value={markerSize}/>
                <TextField step="0.1" type="number" min="0.01" max="1" onKeyPress={this.onMarkerOpacityKeyPress}
                           onChange={this.onMarkerOpacityChange} label="Marker Opacity"
                           className={classes.formControl} value={markerOpacity}/>

                <FormControlLabel
                    control={
                        <Switch
                            checked={binValues}
                            value={'binPlot'}
                            onChange={this.handleBinValuesChange}
                        />
                    }
                    label="Bin Plot"
                />
                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="color-scheme">Color Scheme</InputLabel>
                    <ColorSchemeSelector/>
                </FormControl>
                {binValues &&
                <TextField max="1000" min="20" step="100" onKeyPress={this.onNumberOfBinsKeyPress}
                           value={numberOfBins}
                           onChange={this.onNumberOfBinsChange} label="# Bins Per Axis"
                           className={classes.formControl}/>}

                {binValues && <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="summary">Bin Summary</InputLabel>
                    <Select
                        className={classes.select}
                        input={<Input id="summary"/>}
                        onChange={this.onBinSummaryChange}
                        value={binSummary}
                    >
                        {summaryOptions.map(c => (
                            <MenuItem key={c.value} value={c.value}>
                                <ListItemText primary={c.label}/>
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>}

            </div>
        );
    }
}

const mapStateToProps = state => {
    return {
        selectedFeatures: state.features,
        selectedGroupBy: state.groupBy,
        viewName: state.viewName,
        view3d: state.view3d,
        numberOfBins: state.numberOfBinsUI,
        markerSize: state.markerSizeUI,
        markerOpacity: state.markerOpacityUI,
        binValues: state.binValues,
        binSummary: state.binSummary,
        features: state.dataset == null ? [] : state.dataset.features,
        views: state.dataset == null ? [] : state.dataset.views,
        obsCat: state.dataset == null ? [] : state.dataset.obsCat,
        obs: state.dataset == null ? [] : state.dataset.obs,
    };
};
const mapDispatchToProps = (dispatch, ownProps) => {
    return {
        handleNumberOfBins: value => {
            dispatch(setNumberOfBins(value));
        },
        handleNumberOfBinsUI: value => {
            dispatch(setNumberOfBinsUI(value));
        },
        handleMarkerSize: value => {
            dispatch(setMarkerSize(value));
        },
        handleMarkerSizeUI: value => {
            dispatch(setMarkerSizeUI(value));
        },
        handleMarkerOpacity: value => {
            dispatch(setMarkerOpacity(value));
        },
        handleMarkerOpacityUI: value => {
            dispatch(setMarkerOpacityUI(value));
        },
        handleBinSummary: value => {
            dispatch(setBinSummary(value));
        },
        handleBinValues: value => {
            dispatch(setBinValues(value));
        },
        handleViewName: value => {
            dispatch(setViewName(value));
        },
        handleFeatures: value => {
            dispatch(setFeatures(value == null ? [] : value));
        },
        handleView3d: value => {
            dispatch(setView3d(value));
        },
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(EmbedForm));



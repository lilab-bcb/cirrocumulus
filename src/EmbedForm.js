import {Switch} from '@material-ui/core';
import Divider from '@material-ui/core/Divider';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import InputLabel from '@material-ui/core/InputLabel';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import withStyles from '@material-ui/core/styles/withStyles';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import DeleteIcon from '@material-ui/icons/Delete';
import React from 'react';
import {connect} from 'react-redux';
import {
    deleteDatasetFilter,
    openDatasetFilter,
    setBinSummary,
    setBinValues,
    setEmbeddingChartSize,
    setFeatures,
    setMarkerOpacity,
    setMarkerOpacityUI,
    setMarkerSize,
    setMarkerSizeUI,
    setNumberOfBins,
    setNumberOfBinsUI,
    setSelectedEmbedding,
    setUnselectedMarkerOpacity,
    setUnselectedMarkerOpacityUI,
    setUnselectedMarkerSize,
    setUnselectedMarkerSizeUI
} from './actions';

import AutocompleteSelect from './AutocompleteSelect';
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

    openDatasetFilter = (filterId) => {
        this.props.handleOpenDatasetFilter(filterId);
    };

    deleteDatasetFilter = (filterId) => {
        this.props.handleDeleteDatasetFilter(filterId);
    };


    onMarkerSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            let markerSize = parseFloat(event.target.value);
            if (markerSize >= 0) {
                this.props.handleMarkerSize(markerSize);
            }
        }
    };

    onEmbeddingChartSizeChange = (event) => {
        this.props.handleEmbeddingChartSize(event.target.value);
    };


    onUnselectedMarkerSizeChange = (event) => {
        this.props.handleUnselectedMarkerSizeUI(event.target.value);
    };

    onUnselectedMarkerSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            let markerSize = parseFloat(event.target.value);
            if (markerSize >= 0) {
                this.props.handleUnselectedMarkerSize(markerSize);
            }
        }
    };

    onMarkerOpacityChange = (event) => {
        this.props.handleMarkerOpacityUI(event.target.value);
    };

    onMarkerOpacityKeyPress = (event) => {
        if (event.key === 'Enter') {
            let opacity = parseFloat(event.target.value);
            if (opacity >= 0 && opacity <= 1) {
                this.props.handleMarkerOpacity(opacity);
            }
        }
    };

    onUnselectedMarkerOpacityChange = (event) => {
        this.props.handleUnselectedMarkerOpacityUI(event.target.value);
    };

    onUnselectedMarkerOpacityKeyPress = (event) => {
        if (event.key === 'Enter') {
            let opacity = parseFloat(event.target.value);
            if (opacity >= 0 && opacity <= 1) {
                this.props.handleUnselectedMarkerOpacity(opacity);
            }
        }
    };

    onNumberOfBinsChange = (event) => {
        this.props.handleNumberOfBinsUI(event.target.value);
    };

    onNumberOfBinsKeyPress = (event) => {
        if (event.key === 'Enter') {
            let numberOfBins = parseInt(event.target.value);
            if (numberOfBins >= 0) {
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
        this.props.handleViewName(event.target.value);
    };


    render() {
        const {classes, datasetFilters, embeddingChartSize, unselectedMarkerSize, selectedFeatures, selectedGroupBy, selectedEmbeddings, numberOfBins, markerSize, markerOpacity, unselectedMarkerOpacity, binValues, binSummary, dataset} = this.props;
        let savedDatasetFilter = this.props.savedDatasetFilter;
        if (savedDatasetFilter == null) {
            savedDatasetFilter = {};
        }
        const features = dataset == null ? [] : dataset.features;
        const availableEmbeddings = dataset == null ? [] : dataset.embeddings;
        const isSummarized = dataset == null ? false : dataset.summary != null;
        const obsCat = dataset == null ? [] : dataset.obsCat;
        const obs = dataset == null ? [] : dataset.obs;
        const summaryOptions = [
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
        }];
        const chartSizes = [{label: 'Small', value: 3}, {label: 'Medium', value: 2}, {label: 'Large', value: 1}];
        return (
            <div className={classes.root}>
                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="view">Embedding</InputLabel>
                    <Select
                        className={classes.select}
                        input={<Input id="view"/>}
                        onChange={this.handleViewChange}
                        value={selectedEmbeddings.length === 0 ? '' : selectedEmbeddings[0]}
                        multiple={false}>
                        {availableEmbeddings.map(view => (
                            <MenuItem key={view.name} value={view.name}>
                                <ListItemText primary={view.name}/>
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>


                <FormControl className={classes.formControl}>
                    <AutocompleteSelect label="Features" options={allOptions}
                                        defaultOptions={defaultOptions} value={featureValue}
                                        onChange={this.props.handleFeatures}
                                        isMulti={true}/>
                </FormControl>

                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="chart_size">Chart Size</InputLabel>
                    <Select
                        className={classes.select}
                        input={<Input id="chart_size"/>}
                        onChange={this.onEmbeddingChartSizeChange}
                        value={embeddingChartSize}
                        multiple={false}>
                        {chartSizes.map(item => (
                            <MenuItem key={item.label} value={item.value}>
                                <ListItemText primary={item.label}/>
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>


                <TextField type="number" step="2" min="0" max="30" onKeyPress={this.onMarkerSizeKeyPress}
                           onChange={this.onMarkerSizeChange} label="Marker Size"
                           className={classes.formControl} value={markerSize}/>
                <TextField step="0.1" type="number" min="0" max="1" onKeyPress={this.onMarkerOpacityKeyPress}
                           onChange={this.onMarkerOpacityChange} label="Marker Opacity"
                           className={classes.formControl} value={markerOpacity}/>

                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="color-scheme">Color Scheme</InputLabel>
                    <ColorSchemeSelector/>
                </FormControl>

                {!isSummarized && <FormControlLabel
                    control={
                        <Switch
                            checked={binValues}
                            value={'binPlot'}
                            onChange={this.handleBinValuesChange}
                        />
                    }
                    label="Bin Plot"
                />}

                {!isSummarized && binValues &&
                <TextField max="1000" min="20" step="100" onKeyPress={this.onNumberOfBinsKeyPress}
                           value={numberOfBins}
                           onChange={this.onNumberOfBinsChange} label="# Bins Per Axis"
                           className={classes.formControl}/>}

                {!isSummarized && binValues && <FormControl className={classes.formControl}>
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


                <Divider/>

                <Typography
                    color="textSecondary"
                    display="block"
                    variant="caption"
                >
                    Unselected Chart Properties
                </Typography>
                <TextField type="number" step="2" min="0" max="30" onKeyPress={this.onUnselectedMarkerSizeKeyPress}
                           onChange={this.onUnselectedMarkerSizeChange} label="Unselected Marker Size"
                           className={classes.formControl} value={unselectedMarkerSize}/>
                <TextField step="0.1" type="number" min="0" max="1"
                           onKeyPress={this.onUnselectedMarkerOpacityKeyPress}
                           onChange={this.onUnselectedMarkerOpacityChange} label="Unselected Marker Opacity"
                           className={classes.formControl} value={unselectedMarkerOpacity}/>

                {datasetFilters.length > 0 && <Divider/>}
                {datasetFilters.length > 0 && <h3 style={{marginLeft: 8}}>Saved Filters</h3>}

                <List dense={true}>
                    {datasetFilters.map(item => (
                        <ListItem key={item.id} data-key={item.id} button
                                  selected={item.id === savedDatasetFilter.id}
                                  onClick={e => this.openDatasetFilter(item.id)}>
                            <ListItemText primary={item.name}/>
                            <ListItemSecondaryAction onClick={e => this.deleteDatasetFilter(item.id)}>
                                <IconButton edge="end" aria-label="delete">
                                    <DeleteIcon/>
                                </IconButton>
                            </ListItemSecondaryAction>
                        </ListItem>
                    ))}
                </List>

            </div>
        );
    }
}

const mapStateToProps = state => {
    return {
        binSummary: state.binSummary,
        binValues: state.binValues,
        dataset: state.dataset,
        datasetFilters: state.datasetFilters,
        embeddingChartSize: state.embeddingChartSize,
        markerOpacity: state.markerOpacityUI,
        markerSize: state.markerSizeUI,
        numberOfBins: state.numberOfBinsUI,
        savedDatasetFilter: state.savedDatasetFilter,
        selectedEmbeddings: state.embeddings,
        selectedFeatures: state.features,
        selectedGroupBy: state.groupBy,
        unselectedMarkerOpacity: state.unselectedMarkerOpacityUI,
        unselectedMarkerSize: state.unselectedMarkerSizeUI
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
        handleUnselectedMarkerSize: value => {
            dispatch(setUnselectedMarkerSize(value));
        },
        handleUnselectedMarkerSizeUI: value => {
            dispatch(setUnselectedMarkerSizeUI(value));
        },
        handleMarkerOpacity: value => {
            dispatch(setMarkerOpacity(value));
        },
        handleMarkerOpacityUI: value => {
            dispatch(setMarkerOpacityUI(value));
        },
        handleUnselectedMarkerOpacity: value => {
            dispatch(setUnselectedMarkerOpacity(value));
        },
        handleEmbeddingChartSize: value => {
            dispatch(setEmbeddingChartSize(value));
        },
        handleUnselectedMarkerOpacityUI: value => {
            dispatch(setUnselectedMarkerOpacityUI(value));
        },
        handleBinSummary: value => {
            dispatch(setBinSummary(value));
        },
        handleBinValues: value => {
            dispatch(setBinValues(value));
        },
        handleViewName: value => {
            dispatch(setSelectedEmbedding(value == null ? [] : [value]));
        },
        handleFeatures: value => {
            dispatch(setFeatures(value == null ? [] : value));
        },
        handleOpenDatasetFilter: value => {
            dispatch(openDatasetFilter(value));
        },
        handleDeleteDatasetFilter: value => {
            dispatch(deleteDatasetFilter(value));
        },

    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(EmbedForm));



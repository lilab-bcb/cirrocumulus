import {Switch} from '@material-ui/core';
import Checkbox from '@material-ui/core/Checkbox';
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
            let value = parseInt(event.target.value);
            if (value >= 0) {
                this.props.handleNumberOfBins(event.target.value);
                let embeddings = this.props.embeddings;
                for (let i = 0; i < embeddings.length; i++) {
                    if (!embeddings[i].precomputed) {
                        embeddings[i] = Object.assign(embeddings[i], {nbins: value, _nbins: value});
                    }
                }
                this.props.handleEmbeddings(embeddings.slice(0));
            }
        }
    };

    onBinSummaryChange = (event) => {
        const value = event.target.value;
        this.props.handleBinSummary(value);
        let embeddings = this.props.embeddings;
        for (let i = 0; i < embeddings.length; i++) {
            if (!embeddings[i].precomputed) {
                embeddings[i] = Object.assign(embeddings[i], {agg: value});
            }
        }
        this.props.handleEmbeddings(embeddings.slice(0));
    };

    handleBinValuesChange = (event) => {
        const value = event.target.checked;
        this.props.handleBinValues(value);
        let embeddings = this.props.embeddings;
        for (let i = 0; i < embeddings.length; i++) {
            if (!embeddings[i].precomputed) {
                embeddings[i] = Object.assign(embeddings[i], {bin: value});
            }
        }
        this.props.handleEmbeddings(embeddings.slice(0));
    };

    handleViewChange = (event) => {
        const dataset = this.props.dataset;
        const embeddings = dataset.embeddings;
        const names = event.target.value;
        let selection = [];
        embeddings.forEach(embedding => {
            if (names.indexOf(embedding.name) !== -1) {
                if (!embedding.precomputed) {
                    embedding = Object.assign(embedding, {
                        bin: this.props.binValues,
                        nbins: this.props.numberOfBins,
                        _nbins: this.props.numberOfBinsUI,
                        agg: this.props.binSummary
                    });
                }
                selection.push(embedding);
            }
        });

        this.props.handleEmbeddings(selection);
    };


    render() {
        const {numberOfBinsUI, binValues, binSummary, embeddings, classes, datasetFilters, embeddingChartSize, unselectedMarkerSize, features, groupBy, markerSize, markerOpacity, unselectedMarkerOpacity, dataset} = this.props;

        let savedDatasetFilter = this.props.savedDatasetFilter;
        if (savedDatasetFilter == null) {
            savedDatasetFilter = {};
        }
        const availableFeatures = dataset == null ? [] : dataset.features;
        const availableEmbeddings = dataset == null ? [] : dataset.embeddings;
        const isSummarized = dataset == null ? false : dataset.precomputed != null;
        const obsCat = dataset == null ? [] : dataset.obsCat;
        const obs = dataset == null ? [] : dataset.obs;
        const summaryOptions = [
            {value: 'max', label: 'Maximum'},
            {value: 'mean', label: 'Mean'},
            {value: 'sum', label: 'Sum'}];
        let featureValue = features.concat(groupBy);
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
            options: availableFeatures.map(item => {
                return {label: item, value: item};
            })
        }];
        let defaultOptions = [{
            label: 'Annotations', options: metadataOptions
        }, {
            label: 'Variables',
            options: [{isDisabled: true, label: 'Type to search', value: ''}]
        }];
        const embeddingNames = embeddings.map(e => e.name);
        const chartSizes = [{label: 'Small', value: 3}, {label: 'Medium', value: 2}, {label: 'Large', value: 1}];
        return (
            <div className={classes.root}>
                <FormControl className={classes.formControl}>
                    <InputLabel id="embedding-label">Embedding</InputLabel>
                    <Select
                        className={classes.select}
                        labelId="embedding-label"
                        multiple
                        value={embeddingNames}
                        onChange={this.handleViewChange}
                        input={<Input/>}
                        renderValue={selected => selected.join(', ')}
                    >
                        {availableEmbeddings.map(embedding => (
                            <MenuItem key={embedding.name} value={embedding.name}>
                                <Checkbox checked={embeddingNames.indexOf(embedding.name) > -1}/>
                                <ListItemText primary={embedding.name}/>
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
                           value={numberOfBinsUI}
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
        dataset: state.dataset,
        binValues: state.binValues,
        binSummary: state.binSummary,
        numberOfBins: state.numberOfBins,
        numberOfBinsUI: state.numberOfBinsUI,
        datasetFilters: state.datasetFilters,
        embeddingChartSize: state.embeddingChartSize,
        markerOpacity: state.markerOpacityUI,
        markerSize: state.markerSizeUI,
        savedDatasetFilter: state.savedDatasetFilter,
        embeddings: state.embeddings,
        features: state.features,
        groupBy: state.groupBy,
        unselectedMarkerOpacity: state.unselectedMarkerOpacityUI,
        unselectedMarkerSize: state.unselectedMarkerSizeUI
    };
};
const mapDispatchToProps = (dispatch, ownProps) => {
    return {
        handleEmbeddings: value => {
            dispatch(setSelectedEmbedding(value));
        },
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
        handleFeatures: value => {
            dispatch(setFeatures(value == null ? [] : value));
        },
        handleOpenDatasetFilter: value => {
            dispatch(openDatasetFilter(value));
        },
        handleDeleteDatasetFilter: value => {
            dispatch(deleteDatasetFilter(value));
        }
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(EmbedForm));



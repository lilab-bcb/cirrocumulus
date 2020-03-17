import {Switch} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import Checkbox from '@material-ui/core/Checkbox';
import Divider from '@material-ui/core/Divider';
import MuiExpansionPanel from '@material-ui/core/ExpansionPanel';
import MuiExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import MuiExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
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
import DeleteIcon from '@material-ui/icons/Delete';
import React from 'react';
import {connect} from 'react-redux';
import {
    deleteDatasetFilter,
    exportDatasetFilters,
    getEmbeddingKey,
    handleBrushFilterUpdated,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated,
    openDatasetFilter,
    setBinSummary,
    setBinValues,
    setFeatures,
    setInterpolator,
    setMarkerOpacity,
    setMarkerOpacityUI,
    setNumberOfBins,
    setNumberOfBinsUI,
    setSelectedEmbedding,
    setUnselectedMarkerOpacity,
    setUnselectedMarkerOpacityUI
} from './actions';

import AutocompleteSelect from './AutocompleteSelect';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ColorSchemeSelector from './ColorSchemeSelector';

const styles = theme => ({
    root: {
        display: 'flex',
        flexWrap: 'wrap',
        width: '100%',
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


const ExpansionPanel = withStyles({
    root: {
        border: '1px solid rgba(0, 0, 0, .125)',
        boxShadow: 'none',
        '&:not(:last-child)': {
            borderBottom: 0,
        },
        '&:before': {
            display: 'none',
        },
        '&$expanded': {
            margin: 0,
        },
    },
    expanded: {},
})(MuiExpansionPanel);

const ExpansionPanelSummary = withStyles({
    root: {
        backgroundColor: 'rgba(0, 0, 0, .03)',
        borderBottom: '1px solid rgba(0, 0, 0, .125)',
        marginBottom: -1,
        minHeight: 43,
        '&$expanded': {
            minHeight: 43,
        },
    },
    content: {
        '&$expanded': {
            margin: 0,
        },
    },
    expanded: {},
})(MuiExpansionPanelSummary);

const ExpansionPanelDetails = withStyles(theme => ({
    root: {
        padding: 0,
    },
}))(MuiExpansionPanelDetails);


class EmbedForm extends React.PureComponent {


    openDatasetFilter = (filterId) => {
        this.props.handleOpenDatasetFilter(filterId);
    };

    deleteDatasetFilter = (filterId) => {
        this.props.handleDeleteDatasetFilter(filterId);
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


    handleReverseColors = (event) => {
        const value = event.target.checked;
        const interpolator = this.props.interpolator;
        interpolator.reversed = value;
        this.props.handleInterpolator(Object.assign({}, interpolator));
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

    handleEmbeddingsChange = (event) => {

        const embeddings = event.target.value;
        const selection = [];
        embeddings.forEach(embedding => {
            if (!embedding.precomputed) {
                embedding = Object.assign(embedding, {
                    bin: this.props.binValues,
                    nbins: this.props.numberOfBins,
                    _nbins: this.props.numberOfBinsUI,
                    agg: this.props.binSummary
                });
            }
            selection.push(embedding);
        });
        this.props.handleEmbeddings(selection);
    };


    render() {
        const {
            numberOfBinsUI, interpolator, binValues, binSummary, embeddings, classes, datasetFilters, embeddingData,
            features, groupBy, markerOpacity, datasetFilter,
            featureSummary, shape, nObsSelected, globalFeatureSummary, unselectedMarkerOpacity, dataset,
            handleColorChange, handleMeasureFilterUpdated, handleDimensionFilterUpdated
        } = this.props;

        // for filters we only need one embedding trace per feature
        const traceNames = new Set();
        const filterTraces = [];
        embeddingData.forEach(trace => {
            if (trace.active && trace.name !== '__count' && !traceNames.has(trace.name)) {
                traceNames.add(trace.name);
                filterTraces.push(trace);
            }
        });
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
        const embeddingKeys = embeddings.map(e => getEmbeddingKey(e));
        return (
            <div className={classes.root}>
                <FormControl className={classes.formControl}>
                    <InputLabel id="embedding-label">Embeddings</InputLabel>
                    <Select
                        className={classes.select}
                        labelId="embedding-label"
                        multiple
                        value={embeddings}
                        onChange={this.handleEmbeddingsChange}
                        input={<Input/>}
                        renderValue={selected => selected.map(e => e.name + (e.dimensions === 3 ? ' 3d' : '')).join(', ')}
                    >
                        {availableEmbeddings.map(embedding => (
                            <MenuItem key={getEmbeddingKey(embedding)}
                                      value={embedding}>
                                <Checkbox checked={embeddingKeys.indexOf(getEmbeddingKey(embedding)) !== -1}/>
                                <ListItemText primary={embedding.name + (embedding.dimensions === 3 ? ' 3d' : '')}/>
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>
                <FormControl className={classes.formControl} style={{zIndex: 2000}}>
                    <AutocompleteSelect label="Features" options={allOptions}
                                        defaultOptions={defaultOptions} value={featureValue}
                                        onChange={this.props.handleFeatures}
                                        isMulti={true}/>
                </FormControl>

                <ExpansionPanel defaultExpanded>
                    <ExpansionPanelSummary
                        aria-controls="summary-content"
                        id="summary-header"
                    >
                        <div>Filters</div>
                    </ExpansionPanelSummary>
                    <ExpansionPanelDetails>
                        <div style={{marginLeft: 10}}>
                            {filterTraces.map(traceInfo =>
                                traceInfo.continuous ?
                                    <ColorSchemeLegendWrapper
                                        key={traceInfo.name}
                                        width={140}
                                        showColorScheme={false}
                                        height={30}
                                        handleUpdate={handleMeasureFilterUpdated}
                                        datasetFilter={datasetFilter}
                                        scale={traceInfo.colorScale}
                                        featureSummary={featureSummary}
                                        globalFeatureSummary={globalFeatureSummary}
                                        nObs={shape[0]}
                                        nObsSelected={nObsSelected}
                                        maxHeight={null}
                                        name={traceInfo.name}
                                    /> :
                                    <CategoricalLegend
                                        key={traceInfo.name}
                                        datasetFilter={datasetFilter}
                                        handleClick={handleDimensionFilterUpdated}
                                        handleColorChange={handleColorChange}
                                        name={traceInfo.name}
                                        scale={traceInfo.colorScale}
                                        maxHeight={300}
                                        clickEnabled={true}
                                        nObs={shape[0]}
                                        nObsSelected={nObsSelected}
                                        globalFeatureSummary={globalFeatureSummary}
                                        featureSummary={featureSummary}/>
                            )}

                        </div>
                    </ExpansionPanelDetails>
                </ExpansionPanel>

                <ExpansionPanel defaultExpanded>
                    <ExpansionPanelSummary
                        aria-controls="chart-options-content"
                        id="chart-options-header"
                    >
                        <div>Chart Options</div>
                    </ExpansionPanelSummary>
                    <ExpansionPanelDetails>
                        <div>
                            {/*<FormControl className={classes.formControl}>*/}
                            {/*    <InputLabel htmlFor="chart_size">Chart Size</InputLabel>*/}
                            {/*    <Select*/}
                            {/*        className={classes.select}*/}
                            {/*        input={<Input id="chart_size"/>}*/}
                            {/*        onChange={this.onEmbeddingChartSizeChange}*/}
                            {/*        value={embeddingChartSize}*/}
                            {/*        multiple={false}>*/}
                            {/*        {chartSizes.map(item => (*/}
                            {/*            <MenuItem key={item.label} value={item.value}>*/}
                            {/*                <ListItemText primary={item.label}/>*/}
                            {/*            </MenuItem>*/}
                            {/*        ))}*/}
                            {/*    </Select>*/}
                            {/*</FormControl>*/}


                            {/*<TextField type="text" onKeyPress={this.onMarkerSizeKeyPress}*/}
                            {/*           onChange={this.onMarkerSizeChange} label="Marker Size"*/}
                            {/*           className={classes.formControl} value={markerSize}/>*/}
                            {/*<TextField type="text" onKeyPress={this.onUnselectedMarkerSizeKeyPress}*/}
                            {/*           onChange={this.onUnselectedMarkerSizeChange} label="Unselected Marker Size"*/}
                            {/*           className={classes.formControl} value={unselectedMarkerSize}/>*/}
                            <TextField type="text" onKeyPress={this.onMarkerOpacityKeyPress}
                                       onChange={this.onMarkerOpacityChange} label="Marker Opacity"
                                       className={classes.formControl} value={markerOpacity}/>
                            <TextField type="text"
                                       onKeyPress={this.onUnselectedMarkerOpacityKeyPress}
                                       onChange={this.onUnselectedMarkerOpacityChange} label="Unselected Marker Opacity"
                                       className={classes.formControl} value={unselectedMarkerOpacity}/>
                            <FormControl className={classes.formControl}>
                                <InputLabel htmlFor="color-scheme">Color Scheme</InputLabel>
                                <ColorSchemeSelector/>
                            </FormControl>
                            {<FormControlLabel
                                control={
                                    <Switch
                                        checked={interpolator.reversed || false}
                                        value={'reverseColors'}
                                        onChange={this.handleReverseColors}
                                    />
                                }
                                label="Reverse Colors"
                            />}

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

                            {/*<Typography*/}
                            {/*    color="textSecondary"*/}
                            {/*    display="block"*/}
                            {/*    variant="caption"*/}
                            {/*>*/}
                            {/*    Unselected Chart Properties*/}
                            {/*</Typography>*/}


                        </div>

                    </ExpansionPanelDetails>
                </ExpansionPanel>
                <ExpansionPanel defaultExpanded>
                    <ExpansionPanelSummary
                        aria-controls="filter-options-content"
                        id="filter-options-header"
                    >
                        <div>Saved Filters</div>
                    </ExpansionPanelSummary>
                    <ExpansionPanelDetails>
                        <div>
                            {datasetFilters.length === 0 &&
                            <Box color="text.secondary">No saved filters</Box>}
                            {datasetFilters.length > 0 && <div><List dense={true}>
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
                                <Button onClick={this.props.handleExportDatasetFilters}>Export Filters</Button>
                            </div>}
                        </div>
                    </ExpansionPanelDetails>
                </ExpansionPanel>
            </div>
        );
    }
}

const mapStateToProps = state => {
    return {


        datasetFilter: state.datasetFilter,
        featureSummary: state.featureSummary,
        shape: state.dataset.shape,
        nObsSelected: state.selection.count,
        globalFeatureSummary: state.globalFeatureSummary,
        dataset: state.dataset,
        binValues: state.binValues,
        binSummary: state.binSummary,
        numberOfBins: state.numberOfBins,
        numberOfBinsUI: state.numberOfBinsUI,
        datasetFilters: state.datasetFilters,
        embeddingData: state.embeddingData,
        embeddingChartSize: state.embeddingChartSize,
        interpolator: state.interpolator,
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
        handleDimensionFilterUpdated: (e) => {
            dispatch(handleDimensionFilterUpdated(e));
        },
        handleColorChange: (e) => {
            dispatch(handleColorChange(e));
        },
        handleMeasureFilterUpdated: (e) => {
            dispatch(handleMeasureFilterUpdated(e));
        },
        onSelect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        },
        onDeselect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        },

        handleEmbeddings: value => {
            dispatch(setSelectedEmbedding(value));
        },
        handleInterpolator: value => {
            dispatch(setInterpolator(value));
        },
        handleNumberOfBins: value => {
            dispatch(setNumberOfBins(value));
        },
        handleNumberOfBinsUI: value => {
            dispatch(setNumberOfBinsUI(value));
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
        },
        handleExportDatasetFilters: () => {
            dispatch(exportDatasetFilters());
        }
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(EmbedForm));



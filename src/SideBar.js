import {DialogActions, DialogContentText, InputLabel, Switch, Typography} from '@material-ui/core';

import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Chip from '@material-ui/core/Chip';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Divider from '@material-ui/core/Divider';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormHelperText from '@material-ui/core/FormHelperText';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import Slider from '@material-ui/core/Slider';
import withStyles from '@material-ui/core/styles/withStyles';
import TextField from '@material-ui/core/TextField';
import Tooltip from '@material-ui/core/Tooltip';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import CloudDownloadIcon from '@material-ui/icons/CloudDownload';
import CompareIcon from '@material-ui/icons/Compare';
import DeleteIcon from '@material-ui/icons/Delete';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import FontDownloadRoundedIcon from '@material-ui/icons/FontDownloadRounded';
import HighlightOffIcon from '@material-ui/icons/HighlightOff';
import SaveIcon from '@material-ui/icons/Save';
import {debounce, find, findIndex} from 'lodash';
import memoize from "memoize-one";
import React from 'react';
import {connect} from 'react-redux';
import {AccordionDetailsStyled, AccordionStyled, AccordionSummaryStyled} from './accordion';
import {
    datasetFilterToJson,
    deleteDatasetFilter,
    deleteFeatureSet,
    downloadSelectedIds,
    exportDatasetFilters,
    getDatasetFilterNames,
    getEmbeddingKey,
    getTraceKey,
    handleDomainChange,
    MORE_OPTIONS_DIALOG,
    openDatasetFilter,
    removeDatasetFilter,
    SAVE_DATASET_FILTER_DIALOG,
    SAVE_FEATURE_SET_DIALOG,
    setActiveFeature,
    setChartOptions,
    setChartSize,
    setCombineDatasetFilters,
    setDialog,
    setDistributionPlotOptions,
    setInterpolator,
    setMarkerOpacity,
    setPointSize,
    setSearchTokens,
    setSelectedEmbedding,
    setTab,
    setUnselectedMarkerOpacity,
    submitJob,
    toggleEmbeddingLabel
} from './actions';
import AutocompleteVirtualized from './AutocompleteVirtualized';
import {EditableColorScheme} from './EditableColorScheme';
import {intFormat} from './formatters';
import JobResultOptions from './JobResultOptions';
import {FEATURE_TYPE, getFeatureSets, NATSORT, splitSearchTokens, TRACE_TYPE_META_IMAGE,} from './util';

const pointSizeOptions = [{value: 0.1, label: '10%'}, {value: 0.25, label: '25%'}, {value: 0.5, label: '50%'}, {
    value: 0.75,
    label: '75%'
}, {value: 1, label: '100%'}, {value: 1.5, label: '150%'}, {value: 2, label: '200%'}, {
    value: 3,
    label: '300%'
}, {value: 4, label: '400%'}];
const gallerySizeOptions = [{value: 200, label: 'Extra Small'}, {value: 300, label: 'Small'}, {
    value: 500,
    label: 'Medium'
}, {
    value: 800,
    label: 'Large'
}];

const styles = theme => ({
    root: {
        display: 'flex',
        width: '100%',
        flexDirection: 'column'
    },
    formControl: {
        display: 'block',
        minWidth: 200,
        margin: theme.spacing(0, 1)
    },
    margin: {
        margin: theme.spacing(1, 1)
    },
    select: {
        minWidth: 200,
    },
    toolbar: {
        '& hr': {
            margin: theme.spacing(0, 0.5),
        }
    },
});


const getAnnotationOptions = memoize(
    (obs, obsCat) => {
        const options = [];
        obs.forEach(item => {
            options.push({group: 'Continuous', text: item, id: item});
        });
        obsCat.forEach(item => {
            options.push({group: 'Categorical', text: item, id: item});
        });
        options.sort((item1, item2) => {
            const c = NATSORT(item1.group, item2.group);
            if (c !== 0) {
                return c;
            }
            return NATSORT(item1.text, item2.text);
        });
        return options;
    }
);
const getEmbeddingOptions = memoize(
    (embeddings) => {
        const options = [];
        embeddings.forEach(item => {
            options.push({text: item.name + (item.dimensions === 3 ? ' 3d' : ''), id: getEmbeddingKey(item)});
        });

        options.sort((item1, item2) => {
            return NATSORT(item1.text.toLowerCase(), item2.text.toLowerCase());
        });
        return options;
    }
);
const getMetafeaturesOptions = memoize((items) => {
        if (items) {
            const options = items.slice();
            options.sort(NATSORT);
            return options;
        }
    }
);
const getFeatureSetOptions = memoize((items, categoricalNames) => {
        const options = items.map(item => ({group: item.category, text: item.name, id: item.id}));

        options.forEach(item => {
            let group = item.group;
            let index = group.lastIndexOf(' (');
            if (index !== -1) {
                group = group.substring(0, index);
            }
            let map = categoricalNames[group] || {};
            let newName = map[item.text];
            if (newName != null) {
                item.text = newName;
            }
        });
        options.sort((item1, item2) => {
            const c = NATSORT(item1.group, item2.group);
            if (c !== 0) {
                return c;
            }
            return NATSORT(item1.text, item2.text);
        });
        return options;
    }
);


class SideBar extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            featureSetAnchorEl: null,
            featureSet: null,
            featureSetView: false,
            opacity: props.markerOpacity,
            unselectedOpacity: props.unselectedMarkerOpacity,
            minColor: '',
            maxColor: ''
        };
        this.updateMarkerOpacity = debounce(this.updateMarkerOpacity, 500);
        this.updateUnselectedMarkerOpacity = debounce(this.updateUnselectedMarkerOpacity, 500);
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.markerOpacity !== this.props.markerOpacity) {
            this.setState({opacity: this.props.markerOpacity});
        }
        if (prevProps.unselectedMarkerOpacity !== this.props.unselectedMarkerOpacity) {
            this.setState({unselectedOpacity: this.props.unselectedMarkerOpacity});
        }
        if (prevProps.dataset !== this.props.dataset) {
            this.setState({group1: null, group2: null, group1Count: null, group2Count: null});
        }
        if (prevProps.activeFeature !== this.props.activeFeature) {
            const summary = this.props.activeFeature == null ? null : this.props.globalFeatureSummary[this.props.activeFeature.name];
            if (this.props.activeFeature == null) {
                this.setState({
                    minColor: '',
                    maxColor: ''
                });
            } else {
                const trace = find(this.props.embeddingData, traceInfo => getTraceKey(traceInfo) === this.props.activeFeature.embeddingKey);
                if (trace.type !== TRACE_TYPE_META_IMAGE) {
                    this.setState({
                        minColor: summary.customMin == null ? '' : summary.customMin,
                        maxColor: summary.customMax == null ? '' : summary.customMax
                    });
                } else {
                    this.setState({
                        minColor: summary.customZMin == null ? '' : summary.customZMin,
                        maxColor: summary.customZMax == null ? '' : summary.customZMax
                    });
                }
            }


        }
    }

    onMinUIChange = (value) => {
        this.setState({minColor: value});
    };

    onMaxUIChange = (value) => {
        this.setState({maxColor: value});
    };

    onMinChange = (value) => {
        const {activeFeature, embeddingData, globalFeatureSummary} = this.props;
        const summary = globalFeatureSummary[activeFeature.name];
        const trace = find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);
        if (trace.type !== TRACE_TYPE_META_IMAGE) {
            summary.customMin = isNaN(value) ? undefined : value;
        } else {
            summary.customZMin = isNaN(value) ? undefined : value;
        }
        this.props.onDomain({
            name: activeFeature.name,
            summary: summary
        });
    };

    onMaxChange = (value) => {
        const {activeFeature, embeddingData, globalFeatureSummary} = this.props;
        const summary = globalFeatureSummary[activeFeature.name];
        const trace = find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);
        if (trace.type !== TRACE_TYPE_META_IMAGE) {
            summary.customMax = isNaN(value) ? undefined : value;
        } else {
            summary.customZMax = isNaN(value) ? undefined : value;
        }
        this.props.onDomain({
            name: this.props.activeFeature.name,
            summary: summary
        });
    };

    onMarkerOpacityChange = (event, value) => {
        this.setState({opacity: value});
        this.updateMarkerOpacity(value);
    };

    updateMarkerOpacity = (value) => {
        let opacity = parseFloat(value);
        if (opacity >= 0 && opacity <= 1) {
            this.props.handleMarkerOpacity(opacity);
        }
    };

    onUnselectedMarkerOpacityChange = (event, value) => {
        this.setState({unselectedOpacity: value});
        this.updateUnselectedMarkerOpacity(value);
    };

    updateUnselectedMarkerOpacity = (value) => {
        let opacity = parseFloat(value);
        if (opacity >= 0 && opacity <= 1) {
            this.props.handleUnselectedMarkerOpacity(opacity);
        }
    };

    onReversedChange = (event) => {
        this.props.handleInterpolator(Object.assign({}, this.props.interpolator, {reversed: event.target.checked}));
    };

    openDatasetFilter = (filterId) => {
        this.props.handleOpenDatasetFilter(filterId);
    };

    deleteDatasetFilter = (filterId) => {
        this.props.handleDeleteDatasetFilter(filterId);
    };

    onPointSizeChange = (event) => {
        this.props.handlePointSize(event.target.value);
    };

    onFeaturesChange = (event, value) => {
        this.props.handleSearchTokens(value, FEATURE_TYPE.X);
    };

    onMetafeaturesChange = (event, value) => {
        this.props.handleSearchTokens(value, FEATURE_TYPE.METAFEATURE);
    };

    onObservationsIconClick = (event, option) => {
        this.props.handleEmbeddingLabel(option);
        event.stopPropagation();
    };

    onObservationsChange = (event, value) => {
        let values = [];
        value.forEach(val => {
            if (val.text !== undefined) {
                values.push(val.text);
            } else {
                values.push(val);
            }
        });
        this.props.handleSearchTokens(values, FEATURE_TYPE.OBS);
    };

    onSaveFeatureList = () => {
        this.props.handleDialog(SAVE_FEATURE_SET_DIALOG);
    };

    onEmbeddingsChange = (event, value) => {
        const selection = [];
        const embeddingKeys = this.props.dataset.embeddings.map(item => getEmbeddingKey(item));
        value.forEach(val => {
            const id = val.id !== undefined ? val.id : val;
            const index = embeddingKeys.indexOf(id);
            let embedding = this.props.dataset.embeddings[index];
            if (!embedding.precomputed) {
                embedding = Object.assign(embedding, {
                    bin: this.props.binValues,
                    nbins: this.props.numberOfBins,
                    agg: this.props.binSummary
                });
            }
            selection.push(embedding);
        });
        this.props.handleEmbeddings(selection);
    };

    onFeatureSetsChange = (event, value) => {
        let values = [];
        value.forEach(val => {
            if (val.id !== undefined) {
                values.push(val.id);
            } else {
                values.push(val);
            }
        });
        this.props.handleSearchTokens(values, FEATURE_TYPE.FEATURE_SET);
    };

    onFeatureSetClick = (event, option) => {
        const id = option.id;
        const target = event.target;
        let markers = this.props.markers;

        let featureSet = null;
        for (let i = 0; i < markers.length; i++) {
            if (markers[i].id === id) {
                featureSet = markers[i];
                break;
            }
        }
        if (featureSet == null) {
            console.log(id + ' not found');
        }
        event.stopPropagation();
        this.setState({featureSetAnchorEl: target, featureSet: featureSet});

    };

    onFeatureSetMenuClose = (event) => {
        this.setState({featureSetAnchorEl: null, featureSet: null});
    };

    onViewFeatureSet = (event) => {
        event.stopPropagation();
        this.setState({featureSetAnchorEl: null, featureSetView: true});
    };

    onCloseViewFeatureSetDialog = (event) => {
        event.stopPropagation();
        this.setState({featureSet: null, featureSetView: false});
    };

    onDeleteFeatureSet = (event) => {
        event.stopPropagation();
        let searchTokens = this.props.searchTokens;
        const featureSetId = this.state.featureSet.id;
        let value = searchTokens.filter(token => token.type === FEATURE_TYPE.FEATURE_SET && token.value.id !== featureSetId);
        this.props.handleSearchTokens(value, FEATURE_TYPE.FEATURE_SET);
        this.props.handleDeleteFeatureSet(featureSetId);
        this.setState({featureSetAnchorEl: null, featureSet: null});
    };

    onFilterChipClicked = (event) => {
        this.onFeatureClick(event, event.target.innerText);
    };

    // onNumberOfBinsChange = (event) => {
    //     this.props.handleNumberOfBinsUI(event.target.value);
    //     this.updateNumberOfBins(event.target.value);
    // };
    //
    // updateNumberOfBins = (value) => {
    //
    //     value = parseInt(value);
    //     if (value >= 0) {
    //         this.props.handleNumberOfBins(value);
    //         let embeddings = this.props.embeddings;
    //         for (let i = 0; i < embeddings.length; i++) {
    //             if (!embeddings[i].precomputed) {
    //                 embeddings[i] = Object.assign(embeddings[i], {nbins: value, _nbins: value});
    //             }
    //         }
    //         this.props.handleEmbeddings(embeddings.slice(0));
    //
    //     }
    // };

    // onBinSummaryChange = (event) => {
    //     const value = event.target.value;
    //     this.props.handleBinSummary(value);
    //     let embeddings = this.props.embeddings;
    //     for (let i = 0; i < embeddings.length; i++) {
    //         if (!embeddings[i].precomputed) {
    //             embeddings[i] = Object.assign(embeddings[i], {agg: value});
    //         }
    //     }
    //     this.props.handleEmbeddings(embeddings.slice(0));
    // };
    //
    // handleBinValuesChange = (event) => {
    //     const value = event.target.checked;
    //     this.props.handleBinValues(value);
    //     let embeddings = this.props.embeddings;
    //     for (let i = 0; i < embeddings.length; i++) {
    //         if (!embeddings[i].precomputed) {
    //             embeddings[i] = Object.assign(embeddings[i], {bin: value});
    //         }
    //     }
    //     this.props.handleEmbeddings(embeddings.slice(0));
    // };

    onFeatureClick = (event, option) => {
        event.stopPropagation();
        const value = option.text !== undefined ? option.text : option;
        let galleryTraces = this.props.embeddingData.filter(traceInfo => traceInfo.active);
        for (let i = 0; i < galleryTraces.length; i++) {
            if (galleryTraces[i].name == value) {
                if (this.props.tab !== 'embedding') {
                    this.props.handleTab('embedding');
                }
                this.props.handleActiveFeature({
                    name: galleryTraces[i].name,
                    type: galleryTraces[i].featureType,
                    embeddingKey: getTraceKey(galleryTraces[i])
                });
                break;
            }
        }
    };

    //
    onChartSizeChange = (event) => {
        const value = event.target.value;
        this.props.handleChartSize(value);
    };

    onShowGalleryLabelsChange = (event) => {
        const {chartOptions} = this.props;
        chartOptions.showGalleryLabels = event.target.checked;
        this.props.handleChartOptions(Object.assign({}, chartOptions));
    };

    handleSelectedCellsClick = (event) => {
        event.preventDefault();
        this.props.downloadSelectedIds();
    };

    onJobNameChange = (event) => {
        this.setState({jobName: event.target.value});
    };

    onSubmitJobCancel = () => {
        this.setState({jobName: '', jobParams: null});
    };

    onSubmitJobOK = () => {
        this.state.jobParams.name = this.state.jobName;
        this.props.handleSubmitJob(this.state.jobParams);
        this.setState({jobName: '', jobParams: null});
    };

    onSetGroup = (groupNumber) => {
        const d = {};
        d['group' + groupNumber] = datasetFilterToJson(this.props.dataset, this.props.datasetFilter, this.props.combineDatasetFilters);
        d['group' + groupNumber + 'Count'] = this.props.selection.size;
        this.setState(d);
    };

    onSubmitJob = (jobType) => {
        if (jobType === 'corr') {
            const activeFeature = this.props.activeFeature;
            this.setState({jobName: '', jobParams: {type: 'corr', params: {ref: activeFeature.name}}});
        } else {
            this.setState({
                jobName: '',
                jobParams: {type: 'de', params: {filter: this.state.group1, filter2: this.state.group2}}
            });
        }
    };

    onDatasetFilterChipDeleted = (name) => {
        this.props.removeDatasetFilter(name);
    };

    onDatasetFilterCleared = () => {
        this.props.removeDatasetFilter(null);
    };

    onDatasetFilterSaved = () => {
        this.props.handleDialog(SAVE_DATASET_FILTER_DIALOG);
    };

    handleCombineDatasetFilters = (event) => {
        this.props.handleCombineDatasetFilters(event.target.checked ? 'or' : 'and');
    };

    onChartTypeChange = (event) => {
        this.props.onDistributionPlotOptions({chartType: event.target.value});
    };

    onViolinScaleChange = (event) => {
        this.props.onDistributionPlotOptions({violinScale: event.target.value});
    };

    onViolinShowBoxplot = (event) => {
        this.props.onDistributionPlotOptions({violinShowBoxplot: event.target.checked});
    };


    render() {
        const {
            activeFeature,
            categoricalNames,
            chartOptions,
            chartSize,
            classes,
            combineDatasetFilters,
            dataset,
            datasetFilter,
            datasetFilters,
            distributionPlotOptions,
            embeddingLabels,
            embeddings,
            embeddingData,
            interpolator,
            markers,
            pointSize,
            searchTokens,
            selection,
            serverInfo,
            tab,
            textColor
        } = this.props;
        const primaryTrace = activeFeature == null ? null : find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);
        const chartType = distributionPlotOptions.chartType;
        const datasetFilterKeys = getDatasetFilterNames(datasetFilter);
        datasetFilterKeys.sort(NATSORT);

        let savedDatasetFilter = this.props.savedDatasetFilter;
        if (savedDatasetFilter == null) {
            savedDatasetFilter = {};
        }
        const splitTokens = splitSearchTokens(searchTokens);
        const featureSets = getFeatureSets(markers, splitTokens.featureSets);
        const featureOptions = dataset.features;
        const metafeatureOptions = getMetafeaturesOptions(dataset.metafeatures);
        const obsCat = dataset.obsCat;
        const obs = dataset.obs;
        const annotationOptions = getAnnotationOptions(obs, obsCat);
        const featureSetOptions = getFeatureSetOptions(markers, categoricalNames);
        const embeddingOptions = getEmbeddingOptions(dataset.embeddings);
        const selectedEmbeddings = getEmbeddingOptions(embeddings);
        const dynamic = serverInfo.dynamic;
        const {
            featureSet,
            featureSetAnchorEl,
            featureSetView,
            unselectedOpacity,
            opacity,
            jobParams,
            jobName,
            minColor,
            maxColor
        } = this.state;
        const jobsEnabled = serverInfo.jobs;
        return (
            <div className={classes.root}>
                <Dialog
                    open={jobParams != null}
                    onClose={this.onSubmitJobCancel}
                    aria-labelledby="submit-job-dialog-title"
                    aria-describedby="submit-job-dialog-description"
                >
                    <DialogTitle id="submit-job-dialog-title">Submit Job</DialogTitle>
                    <DialogContent>
                        <DialogContentText id="submit-job-dialog-description">
                            Job Details
                        </DialogContentText>
                        <TextField
                            onChange={this.onJobNameChange}
                            value={jobName}
                            autoFocus
                            margin="dense"
                            label="Job Name"
                            type="text"
                            fullWidth
                        />
                    </DialogContent>
                    <DialogActions>
                        <Button onClick={this.onSubmitJobCancel}>
                            Cancel
                        </Button>
                        <Button disabled={jobName === ''} variant="contained" onClick={e => this.onSubmitJobOK(jobName)}
                                color="primary">
                            Submit
                        </Button>
                    </DialogActions>
                </Dialog>
                <Dialog
                    open={featureSetView}
                    onClose={this.onCloseViewFeatureSetDialog}
                    aria-labelledby="view-set-dialog-title"
                    fullWidth={true}
                    maxWidth={'lg'}
                >
                    <DialogTitle id="view-set-dialog-title">{featureSet ? featureSet.name : ''}</DialogTitle>
                    <DialogContent>
                        <TextField
                            value={featureSet ? featureSet.features.join('\n') : ''}
                            margin="dense"
                            fullWidth
                            readOnly={true}
                            variant="outlined"
                            multiline={true}
                        />
                    </DialogContent>
                </Dialog>
                <Menu
                    id="feature-set-menu"
                    anchorEl={featureSetAnchorEl}
                    open={Boolean(featureSetAnchorEl)}
                    onClose={this.onFeatureSetMenuClose}
                >
                    <MenuItem onClick={this.onViewFeatureSet}>View</MenuItem>
                    <MenuItem divider={true}/>
                    <MenuItem disabled={featureSet && featureSet.readonly}
                              onClick={this.onDeleteFeatureSet}>Delete</MenuItem>


                </Menu>
                <AccordionStyled
                    style={tab === 'embedding' || tab === 'distribution' || tab === 'composition' ? null : {display: 'none'}}
                    defaultExpanded>
                    <AccordionSummaryStyled>
                        <Typography>View</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div>
                            {tab === 'embedding' && embeddingOptions.length > 0 &&
                            <FormControl className={classes.formControl}>

                                <AutocompleteVirtualized label={"Embeddings"}
                                                         options={embeddingOptions}
                                                         getChipTitle={(option) => {
                                                             return option.text;
                                                         }}
                                                         value={selectedEmbeddings}
                                                         getChipText={(option) => option.text}
                                                         getOptionLabel={(option) => option.text}
                                                         getOptionSelected={(option, value) => findIndex(selectedEmbeddings, item => item.id === option.id) !== -1}
                                                         onChange={this.onEmbeddingsChange}
                                />
                            </FormControl>}
                            {featureOptions.length > 0 && <FormControl className={classes.formControl}>
                                <AutocompleteVirtualized onChipClick={this.onFeatureClick}
                                                         label={"Genes/Features"}
                                                         options={featureOptions}
                                                         value={splitTokens.X}
                                                         onChange={this.onFeaturesChange}
                                                         helperText={"Enter or paste list"}
                                />
                            </FormControl>}

                            {annotationOptions.length > 0 && <FormControl className={classes.formControl}>

                                <AutocompleteVirtualized label={"Cell Metadata"}
                                                         options={annotationOptions}
                                                         value={splitTokens.obsCat.concat(splitTokens.obs)}
                                                         onChipClick={this.onFeatureClick}
                                                         groupBy={true}
                                                         getChipIcon={(option) => {
                                                             return splitTokens.obsCat.indexOf(option) !== -1 ?
                                                                 <FontDownloadRoundedIcon
                                                                     onClick={(event) => {
                                                                         this.onObservationsIconClick(event, option);
                                                                     }}
                                                                     title={"Toggle Show/Hide Labels"}
                                                                     style={{
                                                                         marginLeft: 4,
                                                                         marginTop: 0,
                                                                         marginRight: 0,
                                                                         marginBottom: 0
                                                                     }}
                                                                     className={"MuiChip-deleteIcon MuiChip-deleteIconSmall" + (embeddingLabels.indexOf(option) !== -1 ? ' cirro-active' : '')}/> : null;
                                                         }}
                                                         getOptionSelected={(option, value) => option.id === value}
                                                         onChange={this.onObservationsChange}/>
                            </FormControl>}

                            {metafeatureOptions && metafeatureOptions.length > 0 &&
                            <FormControl className={classes.formControl}>
                                <AutocompleteVirtualized onChipClick={this.onMetaFeatureClick}
                                                         label={"Metagenes/Features"}
                                                         options={metafeatureOptions}
                                                         value={splitTokens.metafeatures}
                                                         onChange={this.onMetafeaturesChange}
                                />
                            </FormControl>}
                            {featureSetOptions.length > 0 && <FormControl className={classes.formControl}>

                                <AutocompleteVirtualized label={"Sets"}
                                                         options={featureSetOptions}
                                                         value={featureSets}
                                                         onChipClick={this.onFeatureSetClick}
                                                         getChipTitle={(option) => {
                                                             return option.category + ', ' + option.name;
                                                         }}
                                                         getChipIcon={(option) => {
                                                             return <ArrowDropDownIcon onClick={(event) => {
                                                                 this.onFeatureSetClick(event, option);
                                                             }}/>;
                                                         }}
                                                         groupBy={true}
                                                         onChange={this.onFeatureSetsChange}
                                                         getOptionSelected={(option, value) => option.id === value.id}
                                                         getChipText={option => option.name}/>
                                {dynamic && splitTokens.X.length > 0 &&
                                <Tooltip title={"Save Current Feature List"}>
                                    <IconButton size={'small'} onClick={this.onSaveFeatureList}>
                                        <SaveIcon/>
                                    </IconButton>
                                </Tooltip>
                                }
                            </FormControl>}
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>

                <AccordionStyled style={tab === 'embedding' || tab === 'distribution' ? null : {display: 'none'}}
                                 defaultExpanded>
                    <AccordionSummaryStyled>
                        <Typography>Filters</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div style={{marginLeft: 10, maxHeight: 500}}>
                            <Grid component="label" alignContent={"flex-start"} container alignItems="center"
                                  spacing={0}>
                                <Grid item><InputLabel shrink={true} variant={"standard"}>Combine
                                    Filters</InputLabel></Grid>
                                <Grid item>AND</Grid>
                                <Grid item>
                                    <Switch
                                        size="small"
                                        checked={combineDatasetFilters === 'or'}
                                        onChange={this.handleCombineDatasetFilters}
                                    />
                                </Grid>
                                <Grid item>OR</Grid>
                            </Grid>

                            {datasetFilterKeys.length > 0 && selection.size > 0 &&
                            <>
                                <div style={{marginBottom: 2}}>
                                    {intFormat(selection.size) + " / " + intFormat(dataset.shape[0]) + ": "}
                                    {datasetFilterKeys.map(key => {
                                        return <Chip
                                            onDelete={() => {
                                                this.onDatasetFilterChipDeleted(key);
                                            }}
                                            onClick={this.onFilterChipClicked} size={"small"} variant={"default"}
                                            style={{marginRight: 2, verticalAlign: 'bottom'}}
                                            key={key}
                                            label={key}

                                        />;
                                    })}
                                    <Divider/>
                                    <Grid container alignItems="center" className={classes.toolbar}>
                                        <Tooltip title={"Clear All"}>
                                            <IconButton size={'small'} disabled={datasetFilterKeys.length === 0}
                                                        onClick={this.onDatasetFilterCleared}><HighlightOffIcon/></IconButton>
                                        </Tooltip>
                                        {dynamic && <Tooltip title={"Save Filter"}>
                                            <IconButton size={'small'} disabled={datasetFilterKeys.length === 0}
                                                        onClick={this.onDatasetFilterSaved}><SaveIcon/></IconButton>
                                        </Tooltip>}
                                        <Tooltip title={"Download Selected IDs"}>
                                            <IconButton size={'small'} disabled={datasetFilterKeys.length === 0}
                                                        onClick={this.handleSelectedCellsClick}><CloudDownloadIcon/></IconButton>
                                        </Tooltip>
                                    </Grid>
                                </div>
                            </>
                            }
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>

                {jobsEnabled && <AccordionStyled style={tab === 'embedding' ? null : {display: 'none'}}
                                                 defaultExpanded>
                    <AccordionSummaryStyled>
                        <Typography>Analysis</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div className={this.props.classes.margin}>
                            <Tooltip
                                title={"Find differentially expressed features between two groups of cells"}><Typography>Differential
                                Expression</Typography></Tooltip>
                            <ButtonGroup variant="outlined">
                                <Tooltip
                                    title={'Group one' + (this.state.group1Count ? ' (' + intFormat(this.state.group1Count) + ' cells)' : '')}>
                                    <Button size={"small"} disabled={selection.size === 0}
                                            onClick={event => this.onSetGroup(1)}>1</Button>
                                </Tooltip>
                                <Tooltip
                                    title={'Group two' + (this.state.group2Count ? ' (' + intFormat(this.state.group2Count) + ' cells)' : '')}>
                                    <Button size={"small"} disabled={selection.size === 0}
                                            onClick={event => this.onSetGroup(2)}>2</Button>
                                </Tooltip>
                                <Button startIcon={<CompareIcon/>} size={"small"} variant="outlined"
                                        disabled={selection.size === 0 || this.state.group1 == null || this.state.group2 == null}
                                        onClick={event => this.onSubmitJob('de')}>Go</Button>
                            </ButtonGroup>

                            {/*<Tooltip*/}
                            {/*    title={"Find correlated features in selected cells"}><Typography>Correlation</Typography></Tooltip>*/}
                            {/*<Button style={{minWidth: 40}}*/}
                            {/*        disabled={selection.size === 0 || activeFeature.type !== FEATURE_TYPE.X}*/}
                            {/*        size={"small"} variant="outlined"*/}
                            {/*        onClick={event => this.onSubmitJob('corr')}>Go</Button>*/}
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>}
                <AccordionStyled style={tab === 'distribution' ? null : {display: 'none'}} defaultExpanded>
                    <AccordionSummaryStyled>
                        <Typography>Distribution Options</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div style={{marginTop: 8}}>
                            {chartType === 'violin' && <FormControl className={classes.formControl}>
                                <InputLabel id="violin-scale-label">Scale</InputLabel>
                                <Select
                                    className={classes.select}
                                    labelId="violin-scale-label"
                                    value={distributionPlotOptions.violinScale}
                                    onChange={this.onViolinScaleChange}
                                >
                                    <MenuItem value={'area'}>Area</MenuItem>
                                    <MenuItem value={'width'}>Width</MenuItem>
                                </Select>
                                <FormHelperText>If "area", violins have the same area. If "width", violins have the
                                    same
                                    maximum
                                    width.</FormHelperText>
                            </FormControl>}

                            {chartType === 'violin' && <div><FormControlLabel
                                control={
                                    <Switch
                                        value={"violinShowBoxplot"}
                                        checked={distributionPlotOptions.violinShowBoxplot}
                                        onChange={this.onViolinShowBoxplot}
                                    />
                                }
                                label="Show Box Plot"
                            /></div>}

                            <FormControl className={classes.formControl}>
                                <InputLabel id="dist-chart-type-label">Chart Type</InputLabel>
                                <Select
                                    className={classes.select}
                                    labelId="dist-chart-type-label"
                                    value={chartType}
                                    onChange={this.onChartTypeChange}
                                >
                                    <MenuItem value={'dotplot'}>Dot Plot</MenuItem>
                                    <MenuItem value={'heatmap'}>Heatmap</MenuItem>
                                    <MenuItem value={'violin'}>Violin</MenuItem>
                                </Select>
                            </FormControl>
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>

                <AccordionStyled style={tab === 'embedding' ? null : {display: 'none'}} defaultExpanded>
                    <AccordionSummaryStyled>
                        <Typography>Embedding Options</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div>
                            <InputLabel style={{marginLeft: 8, marginTop: 8}} shrink={true}>Marker
                                Opacity</InputLabel>
                            <Slider
                                min={0.0}
                                max={1}
                                step={0.01}
                                style={{marginLeft: 10, width: '86%'}}
                                valueLabelDisplay="auto"
                                value={opacity}
                                onChange={this.onMarkerOpacityChange} aria-labelledby="continuous-slider"/>

                            <InputLabel style={{marginLeft: 8, marginTop: 8}} shrink={true}>Filtered Marker
                                Opacity</InputLabel>
                            <Slider
                                min={0.0}
                                max={1}
                                step={0.01}
                                style={{marginLeft: 10, width: '86%'}}
                                valueLabelDisplay="auto"
                                value={unselectedOpacity}
                                onChange={this.onUnselectedMarkerOpacityChange}
                                aria-labelledby="continuous-slider"/>


                            <FormControl className={classes.formControl}>
                                <InputLabel htmlFor="point_size">Marker Size</InputLabel>
                                <Select
                                    className={classes.select}
                                    input={<Input id="point_size"/>}
                                    onChange={this.onPointSizeChange}
                                    value={pointSize}
                                    multiple={false}>
                                    {pointSizeOptions.map(item => (
                                        <MenuItem key={item.label} value={item.value}>
                                            <ListItemText primary={item.label}/>
                                        </MenuItem>
                                    ))}
                                </Select>
                            </FormControl>

                            <FormControl className={classes.formControl}>
                                <InputLabel htmlFor="chart_size">Gallery Chart Size</InputLabel>
                                <Select
                                    className={classes.select}
                                    input={<Input id="chart_size"/>}
                                    onChange={this.onChartSizeChange}
                                    value={chartSize}
                                    multiple={false}>
                                    {gallerySizeOptions.map(item => (
                                        <MenuItem key={item.label} value={item.value}>
                                            <ListItemText primary={item.label}/>
                                        </MenuItem>
                                    ))}
                                </Select>
                            </FormControl>

                            <div><FormControlLabel
                                control={
                                    <Switch
                                        checked={chartOptions.showGalleryLabels}
                                        onChange={this.onShowGalleryLabelsChange}
                                    />
                                }
                                label="Gallery Labels"
                            /></div>

                            <div className={this.props.classes.margin}>
                                <EditableColorScheme
                                    textColor={textColor}
                                    interpolator={interpolator}
                                    domain={primaryTrace && primaryTrace.continuous && primaryTrace.name !== '__count' ? primaryTrace.colorScale.domain() : null}
                                    min={minColor}
                                    max={maxColor}
                                    onMinChange={this.onMinChange}
                                    onMaxChange={this.onMaxChange}
                                    onMinUIChange={this.onMinUIChange}
                                    onMaxUIChange={this.onMaxUIChange}
                                    onInterpolator={this.props.handleInterpolator}/>
                            </div>
                            <Button size={'small'}
                                    aria-label="More Options" onClick={this.props.onMoreOptions}>More Options...
                            </Button>

                        </div>

                    </AccordionDetailsStyled>
                </AccordionStyled>
                {dynamic &&
                <AccordionStyled
                    style={tab === 'embedding' || tab === 'distribution' || tab === 'composition' ? null : {display: 'none'}}
                    defaultExpanded>
                    <AccordionSummaryStyled expandIcon={<ExpandMoreIcon/>}>
                        <Typography>Saved Filters</Typography>
                    </AccordionSummaryStyled>
                    <AccordionDetailsStyled>
                        <div style={{width: '100%'}}>
                            {datasetFilters.length === 0 &&
                            <Box color="text.secondary" style={{paddingLeft: 10}}>No saved filters</Box>}
                            {datasetFilters.length > 0 &&
                            <>
                                <List dense={true}>
                                    {datasetFilters.map(item => (
                                        <ListItem key={item.id} data-key={item.id} button
                                                  selected={item.id === savedDatasetFilter.id}
                                                  onClick={e => this.openDatasetFilter(item.id)}>
                                            <ListItemText primary={item.name}/>
                                            <ListItemSecondaryAction
                                                onClick={e => this.deleteDatasetFilter(item.id)}>
                                                <IconButton edge="end" aria-label="delete">
                                                    <DeleteIcon/>
                                                </IconButton>
                                            </ListItemSecondaryAction>
                                        </ListItem>
                                    ))}
                                </List>
                                <div className={this.props.classes.margin}>
                                    <Divider/>
                                    <Tooltip title={"Export Filters"}>
                                        <IconButton size={'small'}
                                                    onClick={this.props.handleExportDatasetFilters}><CloudDownloadIcon/></IconButton>
                                    </Tooltip>

                                </div>
                            </>}
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>}
                <AccordionStyled style={tab === 'results' ? null : {display: 'none'}}
                                 defaultExpanded>
                    <AccordionDetailsStyled>
                        <div className={this.props.classes.margin}>
                            <JobResultOptions/>
                        </div>
                    </AccordionDetailsStyled>
                </AccordionStyled>
            </div>
        );
    }
}

const mapStateToProps = state => {
        return {
            activeFeature: state.activeFeature,
            categoricalNames: state.categoricalNames,
            chartSize: state.chartSize,
            chartOptions: state.chartOptions,
            combineDatasetFilters: state.combineDatasetFilters,
            dataset: state.dataset,
            datasetFilter: state.datasetFilter,
            datasetFilters: state.datasetFilters,
            distributionPlotOptions: state.distributionPlotOptions,
            embeddingData: state.embeddingData,
            embeddingLabels: state.embeddingLabels,
            embeddings: state.embeddings,
            globalFeatureSummary: state.globalFeatureSummary,
            interpolator: state.interpolator,
            markerOpacity: state.markerOpacity,
            markers: state.markers,
            numberOfBins: state.numberOfBins,
            pointSize: state.pointSize,
            savedDatasetFilter: state.savedDatasetFilter,
            searchTokens: state.searchTokens,
            selection: state.selection,
            serverInfo: state.serverInfo,
            tab: state.tab,
            textColor: state.textColor,
            unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
        return {
            handleDialog: (value) => {
                dispatch(setDialog(value));
            },
            handleTab: (value) => {
                dispatch(setTab(value));
            },
            handleActiveFeature: (value) => {
                dispatch(setActiveFeature(value));
            },
            handleInterpolator: value => {
                dispatch(setInterpolator(value));
            },
            handleChartSize: (value) => {
                dispatch(setChartSize(value));
            },
            handleChartOptions: (value) => {
                dispatch(setChartOptions(value));
            },
            onDomain: (value) => {
                dispatch(handleDomainChange(value));
            },
            handleCombineDatasetFilters: (value) => {
                dispatch(setCombineDatasetFilters(value));
            },
            downloadSelectedIds: () => {
                dispatch(downloadSelectedIds());
            },
            removeDatasetFilter: (filter) => {
                dispatch(removeDatasetFilter(filter));
            },
            handleEmbeddings: value => {
                dispatch(setSelectedEmbedding(value));
            },
            handlePointSize: value => {
                dispatch(setPointSize(value));
            },
            handleMarkerOpacity: value => {
                dispatch(setMarkerOpacity(value));
            },
            handleEmbeddingLabel: value => {
                dispatch(toggleEmbeddingLabel(value));
            },
            handleUnselectedMarkerOpacity: value => {
                dispatch(setUnselectedMarkerOpacity(value));
            },
            onMoreOptions: () => {
                dispatch(setDialog(MORE_OPTIONS_DIALOG));
            },
            handleSearchTokens: (value, type) => {
                dispatch(setSearchTokens(value == null ? [] : value, type));
            },
            handleOpenDatasetFilter: value => {
                dispatch(openDatasetFilter(value));
            },
            handleDeleteDatasetFilter: value => {
                dispatch(deleteDatasetFilter(value));
            },
            handleExportDatasetFilters: () => {
                dispatch(exportDatasetFilters());
            },
            handleDeleteFeatureSet: value => {
                dispatch(deleteFeatureSet(value));
            },
            handleSubmitJob: value => {
                dispatch(submitJob(value));
            },
            onDistributionPlotOptions: (payload) => {
                dispatch(setDistributionPlotOptions(payload));
            },

        };
    }
;

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(SideBar));



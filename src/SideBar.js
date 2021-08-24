import {DialogActions, DialogContentText, InputLabel, Switch, Typography} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormHelperText from '@material-ui/core/FormHelperText';
import IconButton from '@material-ui/core/IconButton';
import Input from '@material-ui/core/Input';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import MenuItem from '@material-ui/core/MenuItem';
import Menu from '@material-ui/core/Menu';
import Select from '@material-ui/core/Select';
import Slider from '@material-ui/core/Slider';
import TextField from '@material-ui/core/TextField';
import Tooltip from '@material-ui/core/Tooltip';
import DeleteIcon from '@material-ui/icons/Delete';
import {debounce, find} from 'lodash';
import React from 'react';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import {
    datasetFilterToJson,
    deleteFeatureSet,
    deleteView,
    downloadSelectedIds,
    getTraceKey,
    handleDomainChange,
    openView,
    removeDatasetFilter,
    SAVE_DATASET_FILTER_DIALOG,
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
    setUnselectedMarkerOpacity, setUnselectedPointSize,
    submitJob,
    toggleEmbeddingLabel
} from './actions';
import {EditableColorScheme} from './EditableColorScheme';
import {intFormat} from './formatters';
import JobResultOptions from './JobResultOptions';
import {copyToClipboard, SERVER_CAPABILITY_JOBS, SERVER_CAPABILITY_SAVE_LINKS, TRACE_TYPE_META_IMAGE} from './util';
import ExplorePanel from "./ExplorePanel";
import Link from "@material-ui/core/Link";
import withStyles from "@material-ui/core/styles/withStyles";
import {connect} from 'react-redux';
import Divider from "@material-ui/core/Divider";
import LinkIcon from "@material-ui/icons/Link";

const pointSizeOptions = [{value: 0.1, label: '10%'}, {value: 0.25, label: '25%'}, {value: 0.5, label: '50%'}, {
    value: 0.75,
    label: '75%'
}, {value: 1, label: '100%'}, {value: 1.5, label: '150%'}, {value: 2, label: '200%'}, {
    value: 3,
    label: '300%'
}, {value: 4, label: '400%'}, {value: 5, label: '500%'}];
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
        flexDirection: 'column',
        margin: theme.spacing(0, 0.5)
    },
    title: {textTransform: 'uppercase'},
    formControl: {
        display: 'block',
        minWidth: 200,
        maxWidth: 200
    },
    select: {
        minWidth: 200
    },
    slider: {}
});


class SideBar extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            opacity: props.markerOpacity,
            unselectedOpacity: props.unselectedMarkerOpacity,
            minColor: '',
            maxColor: '',
            labelFontSize: props.chartOptions.labelFontSize,
            labelStrokeWidth: props.chartOptions.labelStrokeWidth
        };

        this.onLabelFontSizeUpdate = debounce(this.onLabelFontSizeUpdate, 500);
        this.onLabelStrokeWidthUpdate = debounce(this.onLabelStrokeWidthUpdate, 500);
        this.updateMarkerOpacity = debounce(this.updateMarkerOpacity, 500);
        this.updateUnselectedMarkerOpacity = debounce(this.updateUnselectedMarkerOpacity, 500);
    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.chartOptions.labelFontSize !== this.props.chartOptions.labelFontSize) {
            this.setState({
                labelFontSize: this.props.chartOptions.labelFontSize
            });
        }
        if (prevProps.chartOptions.labelStrokeWidth !== this.props.chartOptions.labelStrokeWidth) {
            this.setState({
                labelStrokeWidth: this.props.chartOptions.labelStrokeWidth
            });
        }

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
                if (trace == null) {
                    this.setState({
                        minColor: '',
                        maxColor: ''
                    });
                } else if (trace.type !== TRACE_TYPE_META_IMAGE) {
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

    onLabelFontSizeUpdate = (value) => {
        if (!isNaN(value) && value > 0) {
            this.props.chartOptions.labelFontSize = value;
            this.props.handleChartOptions(this.props.chartOptions);
            this.setState({labelFontSize: value});
        }
    };

    onLabelFontSize = (event) => {
        this.setState({labelFontSize: event.target.value});
        this.onLabelFontSizeUpdate(event.target.value);
    };

    onLabelStrokeWidth = (event) => {
        this.setState({labelStrokeWidth: event.target.value});
        this.onLabelStrokeWidthUpdate(event.target.value);
    };

    onLabelStrokeWidthUpdate = (value) => {
        if (!isNaN(value) && value >= 0) {
            this.props.chartOptions.labelStrokeWidth = value;
            this.props.handleChartOptions(this.props.chartOptions);
            this.setState({labelStrokeWidth: value});
        }

    };

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
        value = 0.01;
        console.log(value);
        this.setState({unselectedOpacity: value});
        this.updateUnselectedMarkerOpacity(value);
    };

    updateUnselectedMarkerOpacity = (value) => {
        let opacity = parseFloat(value);
        if (opacity >= 0 && opacity <= 1) {
            this.props.handleUnselectedMarkerOpacity(opacity);
        }
    };

    openView = (id) => {
        this.props.handleOpenView(id);
    };

    deleteView = (id) => {
        this.props.handleDeleteView(id);
    };

    copyView = (id) => {
        let linkText = window.location.protocol + '//' + window.location.host + window.location.pathname;
        linkText += '#q=' + encodeURIComponent(JSON.stringify({link: id}));
        copyToClipboard(linkText);
    };

    onPointSizeChange = (event) => {
        this.props.handlePointSize(event.target.value);
    };

    onUnselectedPointSizeChange = (event) => {
        this.props.handleUnselectedPointSize(event.target.value);
    };


    onChartSizeChange = (event) => {
        const value = event.target.value;
        this.props.handleChartSize(value);
    };

    onShowGalleryLabelsChange = (event) => {
        const {chartOptions} = this.props;
        chartOptions.showGalleryLabels = event.target.checked;
        this.props.handleChartOptions(Object.assign({}, chartOptions));
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

    onSubmitJob = (jobType, version) => {
        if (jobType === 'corr') {
            const activeFeature = this.props.activeFeature;
            this.setState({jobName: '', jobParams: {type: 'corr', params: {ref: activeFeature.name}}});
        } else {
            this.setState({
                jobName: '',
                jobParams: {
                    type: jobType,
                    params: {version: version, filter: this.state.group1, filter2: this.state.group2}
                }
            });
        }
    };


    onViewSaved = () => {
        this.props.handleDialog(SAVE_DATASET_FILTER_DIALOG);
    };

    onDatasetViewSaved = () => {
        this.props.handleDialog(SAVE_DATASET_FILTER_DIALOG);
    };





    render() {
        const {
            activeFeature,
            chartOptions,
            chartSize,
            classes,
            datasetViews,
            embeddingData,
            interpolator,
            pointSize,
            selection,
            serverInfo,
            tab,
            textColor,
            unselectedPointSize
        } = this.props;
        const primaryTrace = activeFeature == null ? null : find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);


        const {
            unselectedOpacity,
            opacity,
            jobParams,
            jobName,
            minColor,
            maxColor
        } = this.state;

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

                <ExplorePanel/>
                {serverInfo.capabilities.has(SERVER_CAPABILITY_JOBS) &&
                <div style={tab === 'embedding' ? null : {display: 'none'}}>
                    <Tooltip
                        title={"Compare two groups of cells"}>
                        <Typography
                            component={"h1"} className={classes.title}>Compare</Typography>
                    </Tooltip>
                    <ButtonGroup variant="outlined" disabled={selection.size === 0}>
                        <Tooltip
                            title={'Group one' + (this.state.group1Count ? ' (' + intFormat(this.state.group1Count) + ' cells)' : '')}>
                            <Button size={"small"}
                                    onClick={event => this.onSetGroup(1)}>1</Button>
                        </Tooltip>
                        <Tooltip
                            title={'Group two' + (this.state.group2Count ? ' (' + intFormat(this.state.group2Count) + ' cells)' : '')}>
                            <Button size={"small"}
                                    onClick={event => this.onSetGroup(2)}>2</Button>
                        </Tooltip>


                    </ButtonGroup>

                    <Button size={"small"} variant={"outlined"}
                            endIcon={<ArrowDropDownIcon/>}
                            disabled={this.state.group1 == null || this.state.group2 == null}
                            onClick={event => this.setState({compareMenu: event.currentTarget})}>
                        GO</Button>
                    <Menu
                        variant={"menu"}
                        id="compare-menu"
                        anchorEl={this.state.compareMenu}
                        open={Boolean(this.state.compareMenu)}
                        onClose={event => this.setState({compareMenu: null})}
                    >
                        {this.props.compareActions.map(action => <MenuItem title={action.title}
                                                                           key={action.jobType}
                                                                           onClick={event => {
                                                                               this.onSubmitJob(action.jobType, action.version);
                                                                               this.setState({compareMenu: null});
                                                                           }}>{action.title}</MenuItem>)}

                    </Menu>
                    {/*<Tooltip*/}
                    {/*    title={"Find correlated features in selected cells"}><Typography>Correlation</Typography></Tooltip>*/}
                    {/*<Button style={{minWidth: 40}}*/}
                    {/*        disabled={selection.size === 0 || activeFeature.type !== FEATURE_TYPE.X}*/}
                    {/*        size={"small"} variant="outlined"*/}
                    {/*        onClick={event => this.onSubmitJob('corr')}>Go</Button>*/}
                </div>}
                {/*<div*/}
                {/*    style={tab === 'distribution' ? null : {display: 'none'}}>*/}
                {/*    <Divider/>*/}
                {/*    <Typography gutterBottom={false} component={"h1"}*/}
                {/*                className={classes.title}>View</Typography>*/}


                {/*</div>*/}

                <div style={tab === 'embedding' ? null : {display: 'none'}}>
                    <Divider/>
                    <Typography gutterBottom={false} component={"h1"} className={classes.title}>View</Typography>
                    <InputLabel style={{marginTop: 8}} shrink={true}>Marker
                        Opacity</InputLabel>
                    <Slider
                        className={classes.slider}
                        min={0.0}
                        max={1}
                        step={0.01}
                        style={{width: 200}}
                        valueLabelDisplay="auto"
                        value={opacity}
                        onChange={this.onMarkerOpacityChange} aria-labelledby="continuous-slider"/>

                    <InputLabel style={{marginTop: 8}} shrink={true}>Filtered Marker
                        Opacity</InputLabel>
                    <Slider
                        className={classes.slider}
                        min={0.0}
                        max={1}
                        step={0.01}
                        style={{width: 200}}
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
                        <InputLabel htmlFor="filtered_point_size">Filtered Marker Size</InputLabel>
                        <Select
                            className={classes.select}
                            input={<Input id="filtered_point_size"/>}
                            onChange={this.onUnselectedPointSizeChange}
                            value={unselectedPointSize}
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

                    <FormControl className={classes.formControl}>
                        <TextField
                            value={this.state.labelFontSize}
                            onChange={this.onLabelFontSize}
                            size="small"
                            InputLabelProps={{shrink: true}}
                            fullWidth
                            label="Label Font Size"
                        />
                    </FormControl>

                    <FormControl className={classes.formControl}>
                        <TextField
                            value={this.state.labelStrokeWidth}
                            onChange={this.onLabelStrokeWidth}
                            size="small"
                            InputLabelProps={{shrink: true}}
                            fullWidth
                            label="Label Shadow Size"
                        />
                    </FormControl>

                </div>

                {serverInfo.capabilities.has(SERVER_CAPABILITY_SAVE_LINKS) &&
                <div
                    style={tab === 'embedding' || tab === 'distribution' || tab === 'composition' ? null : {display: 'none'}}>
                    <Divider/>
                    <Typography gutterBottom={false} component={"h1"}
                                className={classes.title}>Links</Typography>
                    <Divider/>
                    <FormControl className={classes.formControl}>
                        <Tooltip title={"Save Current Visualization State"}><Link
                            style={{
                                float: 'right',
                                fontSize: '0.75rem'
                            }}
                            onClick={this.onViewSaved}>Save</Link></Tooltip></FormControl>

                    {datasetViews.length === 0 &&
                    <Box color="text.secondary">No saved links</Box>}
                    {datasetViews.length > 0 &&
                    <List dense={true} style={{marginTop: 10}}>
                        {datasetViews.map(item => (
                            <ListItem key={item.id} data-key={item.id} button
                                      onClick={e => this.openView(item.id)}>
                                <ListItemText primary={item.name}/>
                                <ListItemSecondaryAction>
                                    <IconButton edge="end" aria-label="copy" onClick={e => this.copyView(item.id)}>
                                        <LinkIcon/>
                                    </IconButton>
                                    <IconButton edge="end" aria-label="delete" onClick={e => this.deleteView(item.id)}>
                                        <DeleteIcon/>
                                    </IconButton>
                                </ListItemSecondaryAction>
                            </ListItem>
                        ))}
                    </List>
                    }

                </div>}

                <div style={tab === 'results' ? null : {display: 'none'}}>
                    <JobResultOptions/>
                </div>
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
            datasetViews: state.datasetViews,
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
            unselectedPointSize: state.unselectedPointSize,
            unselectedMarkerOpacity: state.unselectedMarkerOpacity
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
            removeDatasetFilter: (filter) => {
                dispatch(removeDatasetFilter(filter));
            },
            handleEmbeddings: value => {
                dispatch(setSelectedEmbedding(value));
            },
            handlePointSize: value => {
                dispatch(setPointSize(value));
            },
            handleUnselectedPointSize: value => {
                dispatch(setUnselectedPointSize(value));
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
            handleSearchTokens: (value, type) => {
                dispatch(setSearchTokens(value == null ? [] : value, type));
            },
            handleOpenView: value => {
                dispatch(openView(value));
            },
            handleDeleteView: value => {
                dispatch(deleteView(value));
            },
            handleDeleteFeatureSet: value => {
                dispatch(deleteFeatureSet(value));
            },
            handleSubmitJob: value => {
                dispatch(submitJob(value));
            },
            onDistributionPlotOptions: (payload) => {
                dispatch(setDistributionPlotOptions(payload));
            }
        };
    }
;

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(SideBar));



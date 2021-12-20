import {InputLabel, Switch, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import FormControl from '@mui/material/FormControl';
import FormControlLabel from '@mui/material/FormControlLabel';
import IconButton from '@mui/material/IconButton';
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemSecondaryAction from '@mui/material/ListItemSecondaryAction';
import ListItemText from '@mui/material/ListItemText';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';
import Slider from '@mui/material/Slider';
import InfoIcon from '@mui/icons-material/Info';
import TextField from '@mui/material/TextField';
import Tooltip from '@mui/material/Tooltip';
import DeleteIcon from '@mui/icons-material/Delete';
import {debounce, find} from 'lodash';
import React from 'react';
import {
    deleteView,
    getTraceKey,
    handleDomainChange,
    openView,
    SAVE_DATASET_FILTER_DIALOG,
    setChartOptions,
    setChartSize,
    setDialog,
    setInterpolator,
    setMarkerOpacity,
    setPointSize,
    setUnselectedMarkerOpacity,
    setUnselectedPointSize
} from './actions';
import {EditableColorScheme} from './EditableColorScheme';
import JobResultOptions from './JobResultOptions';
import {copyToClipboard, REACT_MD_OVERRIDES, SERVER_CAPABILITY_LINKS, TRACE_TYPE_META_IMAGE} from './util';
import Link from "@mui/material/Link";
import withStyles from '@mui/styles/withStyles';
import {connect} from 'react-redux';
import Divider from "@mui/material/Divider";
import LinkIcon from "@mui/icons-material/Link";
import ExplorePanel from './ExplorePanel';
import JobPanel from './JobPanel';
import Popover from '@mui/material/Popover';
import ReactMarkdown from 'markdown-to-jsx';

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
        minWidth: 216,
        maxWidth: 216,
        marginBottom: theme.spacing(1)
    },
    select: {
        minWidth: 216
    }
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
            labelStrokeWidth: props.chartOptions.labelStrokeWidth,
            selectedViewEl: null,
            selectedView: null
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

    onInterpolator = (value) => {
        const {activeFeature} = this.props;
        this.props.handleInterpolator({featureType: activeFeature.type, value: value});
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

    viewDetails = (event, id) => {
        event.stopPropagation();
        this.setState({
            selectedViewEl: event.currentTarget,
            selectedView: find(this.props.datasetViews, item => item.id === id)
        });
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
        this.props.handleChartOptions(chartOptions);
    };


    handleCloseViewDetails = () => {
        this.setState({selectedViewEl: null, selectedView: null});
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
            compareActions,
            datasetViews,
            embeddingData,
            interpolator,
            pointSize,
            serverInfo,
            tab,
            unselectedPointSize
        } = this.props;
        const primaryTrace = activeFeature == null ? null : find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);
        const activeInterpolator = activeFeature == null ? null : interpolator[activeFeature.type];

        const {
            unselectedOpacity,
            opacity,
            minColor,
            maxColor,
            selectedViewEl,
            selectedView
        } = this.state;

        return (<div className={classes.root}>
                <ExplorePanel/>
                <JobPanel compareActions={compareActions}/>
                {selectedView && <Popover
                    id={"view-details"}
                    open={Boolean(selectedViewEl)}
                    anchorEl={selectedViewEl}
                    onClose={this.handleCloseViewDetails}
                    anchorOrigin={{
                        vertical: 'bottom',
                        horizontal: 'center'
                    }}
                    transformOrigin={{
                        vertical: 'top',
                        horizontal: 'center'
                    }}
                >
                    <div style={{width: 500}}>
                        {selectedView.last_updated &&
                            <Typography>Last Updated: {new Date(selectedView.last_updated).toDateString()}</Typography>}
                        {selectedView.notes &&
                            <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                           children={selectedView.notes}/>}
                    </div>
                </Popover>}
                <div style={tab === 'embedding' ? null : {display: 'none'}}>
                    <Divider/>
                    <Typography gutterBottom={true} component={"h1"} className={classes.title}>View</Typography>
                    <InputLabel shrink={true}>Opacity</InputLabel>
                    <Slider
                        min={0.0}
                        max={1}
                        step={0.01}
                        sx={{width: 190, marginLeft: 1}}
                        valueLabelDisplay="auto"
                        value={opacity}
                        onChange={this.onMarkerOpacityChange} aria-labelledby="continuous-slider"/>

                    <InputLabel shrink={true}>Filtered Marker Opacity</InputLabel>
                    <Slider
                        min={0.0}
                        max={1}
                        step={0.01}
                        sx={{width: 190, marginLeft: 1}}
                        valueLabelDisplay="auto"
                        value={unselectedOpacity}
                        onChange={this.onUnselectedMarkerOpacityChange}
                        aria-labelledby="continuous-slider"/>


                    <FormControl className={classes.formControl}>
                        <InputLabel htmlFor="point_size">Marker Size</InputLabel>
                        <Select
                            label={"Marker Size"}
                            labelId={"point_size"}
                            className={classes.select}
                            size={"small"}
                            onChange={this.onPointSizeChange}
                            value={pointSize}
                            multiple={false}>
                            {pointSizeOptions.map(item => (
                                <MenuItem key={item.label} value={item.value}>{item.label}</MenuItem>
                            ))}
                        </Select>
                    </FormControl>

                    <FormControl className={classes.formControl}>
                        <InputLabel htmlFor="filtered_point_size">Filtered Marker Size</InputLabel>
                        <Select
                            label={"Filtered Marker Size"}
                            labelId={"filtered_point_size"}
                            size={"small"}
                            className={classes.select}
                            onChange={this.onUnselectedPointSizeChange}
                            value={unselectedPointSize}
                            multiple={false}>
                            {pointSizeOptions.map(item => (
                                <MenuItem key={item.label} value={item.value}>{item.label}</MenuItem>
                            ))}
                        </Select>
                    </FormControl>

                    <FormControl className={classes.formControl}>
                        <InputLabel htmlFor="chart_size">Gallery Chart Size</InputLabel>
                        <Select
                            label={"Gallery Chart Size"}
                            labelId={"chart_size"}
                            size={"small"}
                            className={classes.select}
                            onChange={this.onChartSizeChange}
                            value={chartSize}
                            multiple={false}>
                            {gallerySizeOptions.map(item => (
                                <MenuItem key={item.label} value={item.value}>{item.label}</MenuItem>
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
                        interpolator={activeInterpolator}
                        domain={primaryTrace && primaryTrace.continuous && primaryTrace.name !== '__count' ? primaryTrace.colorScale.domain() : null}
                        min={minColor}
                        max={maxColor}
                        onMinChange={this.onMinChange}
                        onMaxChange={this.onMaxChange}
                        onMinUIChange={this.onMinUIChange}
                        onMaxUIChange={this.onMaxUIChange}
                        onInterpolator={this.onInterpolator}/>

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

                {serverInfo.capabilities.has(SERVER_CAPABILITY_LINKS) &&
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
                                            {(item.notes || item.last_updated) && <IconButton
                                                edge="end"
                                                aria-label="info"
                                                onClick={e => this.viewDetails(e, item.id)}
                                                size="small">
                                                <InfoIcon/>
                                            </IconButton>}
                                            <IconButton
                                                edge="end"
                                                aria-label="copy"
                                                onClick={e => this.copyView(item.id)}
                                                size="small">
                                                <LinkIcon/>
                                            </IconButton>
                                            <IconButton
                                                edge="end"
                                                aria-label="delete"
                                                onClick={e => this.deleteView(item.id)}
                                                size="small">
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
            dataset: state.dataset,
            datasetViews: state.datasetViews,
            embeddingData: state.embeddingData,
            embeddings: state.embeddings,
            globalFeatureSummary: state.globalFeatureSummary,
            interpolator: state.interpolator,
            markerOpacity: state.markerOpacity,
            markers: state.markers,
            pointSize: state.pointSize,
            serverInfo: state.serverInfo,
            tab: state.tab,
            unselectedPointSize: state.unselectedPointSize,
            unselectedMarkerOpacity: state.unselectedMarkerOpacity
        };
    }
;
const mapDispatchToProps = (dispatch) => {
        return {
            handleDialog: (value) => {
                dispatch(setDialog(value));
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
            handlePointSize: value => {
                dispatch(setPointSize(value));
            },
            handleUnselectedPointSize: value => {
                dispatch(setUnselectedPointSize(value));
            },
            handleMarkerOpacity: value => {
                dispatch(setMarkerOpacity(value));
            },

            handleUnselectedMarkerOpacity: value => {
                dispatch(setUnselectedMarkerOpacity(value));
            },
            handleOpenView: value => {
                dispatch(openView(value));
            },
            handleDeleteView: value => {
                dispatch(deleteView(value));
            }


        };
    }
;

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(SideBar));



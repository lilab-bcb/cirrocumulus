import {InputLabel, Switch, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import FormControl from '@mui/material/FormControl';
import FormControlLabel from '@mui/material/FormControlLabel';
import IconButton from '@mui/material/IconButton';
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemText from '@mui/material/ListItemText';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';
import Slider from '@mui/material/Slider';
import TextField from '@mui/material/TextField';
import Tooltip from '@mui/material/Tooltip';
import {debounce, find} from 'lodash';
import React, {useEffect, useMemo, useState} from 'react';
import {
    deleteLink,
    getTraceKey,
    handleDomainChange,
    openLink,
    SAVE_DATASET_FILTER_DIALOG,
    setChartOptions,
    setChartSize,
    setDialog,
    setInterpolator,
    setMarkerOpacity,
    setMessage,
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
import ExplorePanel from './ExplorePanel';
import JobPanel from './JobPanel';
import Popover from '@mui/material/Popover';
import ReactMarkdown from 'markdown-to-jsx';
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import Menu from '@mui/material/Menu';

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


function SideBar(props) {

    const [opacity, setOpacity] = useState(props.markerOpacity);
    const [unselectedOpacity, setUnselectedOpacity] = useState(props.unselectedMarkerOpacity);
    const [labelFontSize, setLabelFontSize] = useState(props.chartOptions.labelFontSize);
    const [labelStrokeWidth, setLabelStrokeWidth] = useState(props.chartOptions.labelStrokeWidth);

    const [minColor, setMinColor] = useState('');
    const [maxColor, setMaxColor] = useState('');
    const [selectedViewEl, ssetSelectedViewEl] = useState(null);
    const [selectedView, setSelectedView] = useState(null);
    const [contextMenu, setContextMenu] = useState(null);
    const [selectedLink, setSelectedLink] = useState(null);

    const {
        activeFeature,
        chartOptions,
        chartSize,
        classes,
        compareActions,
        datasetViews,
        embeddingData,
        globalFeatureSummary,
        interpolator,
        markerOpacity,
        pointSize,
        serverInfo,
        tab,
        unselectedPointSize,
        unselectedMarkerOpacity,
        handleDialog,
        handleInterpolator,
        handleChartSize,
        handleChartOptions,
        onDomain,
        handlePointSize,
        handleUnselectedPointSize,
        handleMarkerOpacity,
        handleUnselectedMarkerOpacity,
        handleOpenView,
        handleDeleteView,
        handleMessage
    } = props;

    const primaryTrace = activeFeature == null ? null : find(embeddingData, trace => getTraceKey(trace) === activeFeature.embeddingKey);
    const activeInterpolator = activeFeature == null ? null : interpolator[activeFeature.type];
    const onLabelFontSizeUpdateDebouncedFunc = useMemo(() => debounce(onLabelFontSizeUpdate, 500), []);
    const onLabelStrokeWidthUpdateDebouncedFunc = useMemo(() => debounce(onLabelStrokeWidthUpdate, 500), []);
    const updateMarkerOpacityDebouncedFunc = useMemo(() => debounce(updateMarkerOpacity, 500), []);
    const updateUnselectedMarkerOpacityDebouncedFunc = useMemo(() => debounce(updateUnselectedMarkerOpacity, 500), []);

    useEffect(() => {
        setLabelFontSize(chartOptions.labelFontSize);
    }, [chartOptions.labelFontSize]);
    useEffect(() => {
        setLabelStrokeWidth(chartOptions.labelStrokeWidth);
    }, [chartOptions.labelStrokeWidth]);
    useEffect(() => {
        setOpacity(markerOpacity);
    }, [markerOpacity]);

    useEffect(() => {
        setUnselectedOpacity(unselectedMarkerOpacity);
    }, [unselectedMarkerOpacity]);

    useEffect(() => {
        const summary = activeFeature == null ? null : globalFeatureSummary[activeFeature.name];
        if (activeFeature == null) {
            setMinColor('');
            setMaxColor('');
        } else {
            const trace = find(embeddingData, trace => getTraceKey(trace) === activeFeature.embeddingKey);
            if (trace == null) {
                setMinColor('');
                setMaxColor('');
            } else if (trace.type !== TRACE_TYPE_META_IMAGE) {
                setMinColor(summary.customMin == null ? '' : summary.customMin);
                setMaxColor(summary.customMax == null ? '' : summary.customMax);
            } else {
                setMinColor(summary.customZMin == null ? '' : summary.customZMin);
                setMaxColor(summary.customZMax == null ? '' : summary.customZMax);
            }
        }
    }, [activeFeature, embeddingData, globalFeatureSummary]);


    function onLinkContextMenu(event, item) {
        event.preventDefault();
        event.stopPropagation();
        setSelectedLink(item);
        setContextMenu({
            mouseX: event.clientX - 2, mouseY: event.clientY - 4
        });

    }

    function onLabelFontSizeUpdate(value) {
        if (!isNaN(value) && value > 0) {
            chartOptions.labelFontSize = value;
            handleChartOptions(chartOptions);
            setLabelFontSize(value);
        }
    }

    function onLabelFontSize(event) {
        setLabelFontSize(event.target.value);
        onLabelFontSizeUpdateDebouncedFunc(event.target.value);
    }

    function onLabelStrokeWidth(event) {
        setLabelStrokeWidth(event.target.value);
        onLabelStrokeWidthUpdateDebouncedFunc(event.target.value);
    }

    function onLabelStrokeWidthUpdate(value) {
        if (!isNaN(value) && value >= 0) {
            chartOptions.labelStrokeWidth = value;
            handleChartOptions(chartOptions);
            setLabelStrokeWidth(value);
        }
    }

    function onInterpolator(value) {
        handleInterpolator({featureType: activeFeature.type, value: value});
    }

    function onMinUIChange(value) {
        setMinColor(value);
    }

    function onMaxUIChange(value) {
        setMaxColor(value);
    }

    function onMinChange(value) {
        const summary = globalFeatureSummary[activeFeature.name];
        const trace = find(embeddingData, trace => getTraceKey(trace) === activeFeature.embeddingKey);
        if (trace.type !== TRACE_TYPE_META_IMAGE) {
            summary.customMin = isNaN(value) ? undefined : value;
        } else {
            summary.customZMin = isNaN(value) ? undefined : value;
        }
        onDomain({
            name: activeFeature.name,
            summary: summary
        });
    }

    function onMaxChange(value) {
        const summary = globalFeatureSummary[activeFeature.name];
        const trace = find(embeddingData, trace => getTraceKey(trace) === activeFeature.embeddingKey);
        if (trace.type !== TRACE_TYPE_META_IMAGE) {
            summary.customMax = isNaN(value) ? undefined : value;
        } else {
            summary.customZMax = isNaN(value) ? undefined : value;
        }
        onDomain({
            name: activeFeature.name,
            summary: summary
        });
    }

    function onMarkerOpacityChange(event, value) {
        setOpacity(value);
        updateMarkerOpacityDebouncedFunc(value);
    }

    function updateMarkerOpacity(value) {
        let opacity = parseFloat(value);
        if (opacity >= 0 && opacity <= 1) {
            handleMarkerOpacity(opacity);
        }
    }

    function onUnselectedMarkerOpacityChange(event, value) {
        setUnselectedOpacity(value);
        updateUnselectedMarkerOpacityDebouncedFunc(value);
    }

    function updateUnselectedMarkerOpacity(value) {
        let opacity = parseFloat(value);
        if (opacity >= 0 && opacity <= 1) {
            handleUnselectedMarkerOpacity(opacity);
        }
    }

    function openView(id) {
        handleOpenView(id);
    }

    function deleteView(id) {
        setContextMenu(null);
        handleDeleteView(id);
    }

    function viewDetails(event, id) {
        event.stopPropagation();
        setContextMenu(null);
        ssetSelectedViewEl(event.currentTarget);
        setSelectedView(find(datasetViews, item => item.id === id));
    }

    function copyView(id) {
        let linkText = window.location.protocol + '//' + window.location.host + window.location.pathname;
        linkText += '#q=' + encodeURIComponent(JSON.stringify({link: id}));
        setContextMenu(null);
        copyToClipboard(linkText);
        handleMessage('Link copied');
    }

    function onPointSizeChange(event) {
        handlePointSize(event.target.value);
    }

    function onUnselectedPointSizeChange(event) {
        handleUnselectedPointSize(event.target.value);
    }

    function onChartSizeChange(event) {
        handleChartSize(event.target.value);
    }

    function onShowGalleryLabelsChange(event) {
        chartOptions.showGalleryLabels = event.target.checked;
        handleChartOptions(chartOptions);
    }

    function handleCloseViewDetails() {
        ssetSelectedViewEl(null);
        setSelectedView(null);
    }

    function onViewSaved() {
        handleDialog(SAVE_DATASET_FILTER_DIALOG);
    }

    return (<div className={classes.root}>
            <Menu
                anchorReference="anchorPosition"
                anchorPosition={contextMenu != null ? {
                    top: contextMenu.mouseY,
                    left: contextMenu.mouseX
                } : undefined}
                open={Boolean(contextMenu)}
                onClose={e => {
                    setContextMenu(null);
                    setSelectedLink(null);
                }}
            >
                {selectedLink && (selectedLink.notes || selectedLink.last_updated) &&
                    <MenuItem onClick={e => viewDetails(e, selectedLink.id)}>Info</MenuItem>}
                {selectedLink && <MenuItem onClick={e => copyView(selectedLink.id)}>Copy</MenuItem>}
                {selectedLink && <MenuItem onClick={e => deleteView(selectedLink.id)}>Delete</MenuItem>}
            </Menu>
            <ExplorePanel/>
            <JobPanel compareActions={compareActions}/>
            {selectedView && <Popover
                id={"view-details"}
                open={Boolean(selectedViewEl)}
                anchorEl={selectedViewEl}
                onClose={handleCloseViewDetails}
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
                    onChange={onMarkerOpacityChange} aria-labelledby="continuous-slider"/>

                <InputLabel shrink={true}>Filtered Marker Opacity</InputLabel>
                <Slider
                    min={0.0}
                    max={1}
                    step={0.01}
                    sx={{width: 190, marginLeft: 1}}
                    valueLabelDisplay="auto"
                    value={unselectedOpacity}
                    onChange={onUnselectedMarkerOpacityChange}
                    aria-labelledby="continuous-slider"/>


                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="point_size">Marker Size</InputLabel>
                    <Select
                        label={"Marker Size"}
                        labelId={"point_size"}
                        className={classes.select}
                        size={"small"}
                        onChange={onPointSizeChange}
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
                        onChange={onUnselectedPointSizeChange}
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
                        onChange={onChartSizeChange}
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
                            onChange={onShowGalleryLabelsChange}
                        />
                    }
                    label="Gallery Labels"
                /></div>


                <EditableColorScheme
                    interpolator={activeInterpolator}
                    domain={primaryTrace && primaryTrace.continuous && primaryTrace.name !== '__count' ? primaryTrace.colorScale.domain() : null}
                    min={minColor}
                    max={maxColor}
                    onMinChange={onMinChange}
                    onMaxChange={onMaxChange}
                    onMinUIChange={onMinUIChange}
                    onMaxUIChange={onMaxUIChange}
                    onInterpolator={onInterpolator}/>

                <FormControl className={classes.formControl}>
                    <TextField
                        value={labelFontSize}
                        onChange={onLabelFontSize}
                        size="small"
                        sx={{width: 130}}
                        InputLabelProps={{shrink: true}}
                        fullWidth
                        label="Label Font Size"
                    />
                </FormControl>

                <FormControl className={classes.formControl}>
                    <TextField
                        value={labelStrokeWidth}
                        onChange={onLabelStrokeWidth}
                        size="small"
                        sx={{width: 130}}
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
                            onClick={onViewSaved}>Save</Link></Tooltip></FormControl>

                    {datasetViews.length === 0 &&
                        <Box color="text.secondary">No saved links</Box>}
                    {datasetViews.length > 0 &&
                        <List dense={true} style={{marginTop: 10}}>
                            {datasetViews.map(item => (
                                <ListItem key={item.id} data-key={item.id} button
                                          onClick={e => openView(item.id)} secondaryAction={<IconButton
                                    edge="end"
                                    disableRipple={true}
                                    onClick={event => onLinkContextMenu(event, item)}
                                    aria-label="menu"
                                    size="small">
                                    <ArrowDropDownIcon/>
                                </IconButton>}>
                                    <ListItemText primary={item.name}/>
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
                dispatch(openLink(value));
            },
            handleDeleteView: value => {
                dispatch(deleteLink(value));
            },
            handleMessage: value => {
                dispatch(setMessage(value));
            }
        };
    }
;

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(SideBar));



import {Divider, IconButton, Menu, Tooltip} from '@material-ui/core';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import Link from '@material-ui/core/Link';
import MenuItem from '@material-ui/core/MenuItem';
import Popover from '@material-ui/core/Popover';
import {withStyles} from '@material-ui/core/styles';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import AccountCircle from '@material-ui/icons/AccountCircle';
import Brightness2Icon from '@material-ui/icons/Brightness3';
import HelpIcon from '@material-ui/icons/Help';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import ReactMarkdown from 'markdown-to-jsx';
import React from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    login,
    logout,
    setChartOptions,
    setDataset,
    setDialog,
    setMessage,
    setSavedDatasetState,
    setTab,
} from './actions';
import {drawerWidth} from './App';
import CirroIcon from './CirroIcon';
import DatasetSelector from './DatasetSelector';
import {intFormat} from './formatters';
import {
    DEFAULT_DARK_MODE,
    DEFAULT_LABEL_FONT_SIZE,
    DEFAULT_LABEL_STROKE_WIDTH,
    DEFAULT_MARKER_OPACITY,
    DEFAULT_SHOW_FOG,
    DEFAULT_UNSELECTED_MARKER_OPACITY
} from "./reducers";
import {FEATURE_TYPE, REACT_MD_OVERRIDES} from './util';


const styles = theme => ({
    root: {
        display: 'flex',
        flexWrap: 'wrap',
        'flex-direction': 'column',
    },
    appBar: {
        width: `calc(100% - ${drawerWidth}px)`,
        marginLeft: drawerWidth,
    },
});
const AntTab = withStyles(theme => ({
    root: {
        minWidth: 50,
        textTransform: 'none',
        fontWeight: theme.typography.fontWeightRegular,
        marginRight: theme.spacing(0),
        '&:hover': {
            color: '#40a9ff',
            opacity: 1,
        },
        '&$selected': {
            color: '#1890ff',
            fontWeight: theme.typography.fontWeightMedium,
        },
        '&:focus': {
            color: '#40a9ff',
        },
    },
    selected: {},
}))(props => <Tab disableRipple {...props} />);

class AppHeader extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            userMenuOpen: false,
            userMenuAnchorEl: null,
            moreMenuOpen: false,
            moreMenuAnchorEl: null,
        };

    }

    handleTabChange = (event, value) => {
        this.props.handleTab(value);
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


    handleUserMenuClose = () => {
        this.setState({userMenuOpen: false});
    };

    handleMoreMenuClose = () => {
        this.setState({moreMenuOpen: false});
    };


    handleHelp = () => {
        this.props.handleDialog(HELP_DIALOG);
    };


    handleUserMenuOpen = (event) => {
        this.setState({userMenuOpen: true, userMenuAnchorEl: event.currentTarget});
    };
    handleMoreMenuOpen = (event) => {
        this.setState({moreMenuOpen: true, moreMenuAnchorEl: event.currentTarget});
    };

    getLinkJson = () => {
        const {
            chartOptions,
            combineDatasetFilters,
            activeFeature,
            dataset,
            embeddingLabels,
            distributionPlotOptions,
            distributionPlotInterpolator,
            embeddings,
            searchTokens,
            datasetFilter,
            interpolator,
            markerOpacity,
            pointSize,
            unselectedMarkerOpacity,
            distributionData
        } = this.props;

        let json = {
            dataset: dataset.id,
            embeddings: embeddings.map(embedding => {
                if (embedding.bin) {
                    embedding = Object.assign({}, embedding);
                    delete embedding._bin;
                    return embedding;
                } else {
                    return {name: embedding.name, dimensions: embedding.dimensions};
                }

            })
        };

        const chartRef = chartOptions.ref;
        if (json.embeddings.length > 0 && chartRef != null && chartRef.scatterPlot) {
            json.camera = chartRef.scatterPlot.getCameraDef();

        }
        if (activeFeature != null) {
            json.activeFeature = activeFeature;
        }
        let jsonChartOptions = {};

        const defaultChartOptions = {
            showFog: DEFAULT_SHOW_FOG, darkMode: DEFAULT_DARK_MODE,
            labelFontSize: DEFAULT_LABEL_FONT_SIZE,
            labelStrokeWidth: DEFAULT_LABEL_STROKE_WIDTH
        };

        for (let key in defaultChartOptions) {
            let value = chartOptions[key];
            if (value !== defaultChartOptions[key]) {
                jsonChartOptions[key] = value;
            }
        }
        if (pointSize !== 1) {
            json.pointSize = pointSize;
        }

        json.chartOptions = jsonChartOptions;


        if (searchTokens.length > 0) {
            json.q = searchTokens;
        }

        let datasetFilterJson = {};
        for (let key in datasetFilter) {
            let filterObject = datasetFilter[key];
            if (Array.isArray(filterObject)) { // brush filter
                const array = [];
                filterObject.forEach(brush => {
                    const brushJson = Object.assign({}, brush);
                    brushJson.indices = Array.from(brush.indices);
                    array.push(brushJson);
                });
                datasetFilterJson[key] = array;
            } else if (filterObject.operation === 'in') {
                datasetFilterJson[key] = {operation: filterObject.operation, value: filterObject.value};
            } else {
                if (filterObject.operation !== '' && !isNaN(filterObject.value) && filterObject.value != null) {
                    datasetFilterJson[key] = {operation: filterObject.operation, value: filterObject.value};
                }
            }
        }
        if (Object.keys(datasetFilterJson).length > 0) {
            json.datasetFilter = datasetFilterJson;
        }
        json.combineDatasetFilters = combineDatasetFilters;
        if (markerOpacity !== DEFAULT_MARKER_OPACITY) {
            json.markerOpacity = markerOpacity;
        }
        if (unselectedMarkerOpacity !== DEFAULT_UNSELECTED_MARKER_OPACITY) {
            json.unselectedMarkerOpacity = unselectedMarkerOpacity;
        }

        if (distributionData && distributionData.length > 0) {
            json.distributionPlotOptions = distributionPlotOptions;
            json.distributionPlotInterpolator = Object.assign({}, distributionPlotInterpolator, {value: null});
        }

        const interpolatorJson = Object.assign({}, interpolator, {value: null});
        json.colorScheme = interpolatorJson;

        if (embeddingLabels.length > 0) {
            json.embeddingLabels = embeddingLabels;
        }
        return json;
    };

    handleDataset = (id) => {
        if (this.props.dataset != null) {
            const savedDatasetState = this.props.savedDatasetState;
            const link = this.getLinkJson();
            link.dataset = null;
            savedDatasetState[this.props.dataset.id] = link;
            this.props.handleSavedDatasetState(savedDatasetState);
        }
        this.props.handleTab('embedding'); // embedding won't render unless visible
        this.props.handleDataset(id);
    };

    copyLink = (event) => {

        let linkText = window.location.protocol + '//' + window.location.host + window.location.pathname;
        linkText += '?q=' + encodeURIComponent(JSON.stringify(this.getLinkJson()));
        const container = document.activeElement;
        const fakeElem = document.createElement('textarea');
        const isRTL = document.documentElement.getAttribute('dir') == 'rtl';
        fakeElem.style.fontSize = '12pt';
        fakeElem.style.border = '0';
        fakeElem.style.padding = '0';
        fakeElem.style.margin = '0';
        // Move element out of screen horizontally
        fakeElem.style.position = 'absolute';
        fakeElem.style[isRTL ? 'right' : 'left'] = '-9999px';
        fakeElem.setAttribute('readonly', '');
        // Move element to the same position vertically
        let yPosition = window.pageYOffset || document.documentElement.scrollTop;
        fakeElem.style.top = yPosition + 'px';
        fakeElem.value = linkText;
        container.appendChild(fakeElem);

        fakeElem.select();
        fakeElem.setSelectionRange(0, fakeElem.value.length);

        document.execCommand('copy');

        const fakeHandlerCallback = (event) => {
            document.activeElement.blur();
            window.getSelection().removeAllRanges();
            container.removeChild(fakeElem);
            container.removeEventListener('click', fakeHandlerCallback);
        };
        this.props.setMessage('Link copied');
        container.addEventListener('click', fakeHandlerCallback);
        this.setState({moreMenuOpen: false});

    };


    onDarkMode = () => {
        this.props.chartOptions.darkMode = !this.props.chartOptions.darkMode;
        this.props.handleChartOptions(this.props.chartOptions);
    };

    handleLogout = () => {
        this.setState({userMenuOpen: false});
        this.props.handleLogout();
    };

    handleImportDataset = (event) => {
        this.props.handleDialog(IMPORT_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };

    handleSettings = (event) => {
        this.props.handleDialog(EDIT_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };

    handleDelete = (event) => {
        this.props.handleDialog(DELETE_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };

    handleShowDatasetDetails = (event) => {
        this.setState({datasetDetailsEl: event.currentTarget});
    };
    handleCloseDatasetDetails = (event) => {
        this.setState({datasetDetailsEl: null});
    };

    render() {
        const {
            dataset,
            distributionData,
            loadingApp,
            jobResults,
            email,
            selection,
            classes,
            searchTokens,
            serverInfo,
            tab,
            user
        } = this.props;


        const datasetDetailsOpen = Boolean(this.state.datasetDetailsEl);
        const shape = dataset != null && dataset.shape != null ? dataset.shape : null;
        const hasSelection = dataset != null && shape != null && shape[0] > 0 && selection.size > 0;
        const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT).map(item => item.value);
        const showNewDataset = user != null && user.importer && !loadingApp.loading;
        const showEditDeleteDataset = dataset !== null && dataset.owner && !loadingApp.loading;
        const showMoreMenu = (showNewDataset || dataset != null) && !loadingApp.loading;
        const isSignedOut = !loadingApp.loading && email == null && serverInfo.clientId !== '';
        return (
            <AppBar position="fixed" color="default" className={classes.appBar}>
                <Toolbar variant="dense" style={{paddingLeft: 6}}>
                    {dataset != null && <Popover
                        id={"dataset-details"}
                        open={datasetDetailsOpen}
                        anchorEl={this.state.datasetDetailsEl}
                        onClose={this.handleCloseDatasetDetails}
                        anchorOrigin={{
                            vertical: 'bottom',
                            horizontal: 'center',
                        }}
                        transformOrigin={{
                            vertical: 'top',
                            horizontal: 'center',
                        }}
                    >

                        <div style={{width: 500, padding: '1em'}}>
                            <Typography className={classes.typography}>
                                {dataset.name}
                            </Typography>
                            <Divider/>
                            {dataset.species && <Typography className={classes.typography}>
                                Species: {dataset.species}
                            </Typography>}

                            {dataset.title && <Typography className={classes.typography}>
                                Title: {dataset.title}
                            </Typography>}
                            {dataset.description &&
                            <>Description: <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                                          children={dataset.description}/></>}
                            {!process.env.REACT_APP_STATIC === 'true' && <Typography className={classes.typography}>
                                URL: {dataset.url}
                            </Typography>}
                        </div>
                    </Popover>
                    }
                    {dataset == null && <Typography variant="h5">
                        <CirroIcon/> Cirrocumulus
                    </Typography>}
                    {dataset && <CirroIcon/>}
                    {dataset &&
                    <Link
                        style={{
                            paddingLeft: 6,
                            whiteSpace: 'nowrap',
                            maxWidth: 300,
                            overflow: 'hidden',
                            textOverflow: 'ellipsis'
                        }}
                        href="#"
                        underline={"none"}
                        onClick={this.handleShowDatasetDetails}
                        aria-owns={this.state.datasetDetailsOpen ? 'dataset-details' : undefined}
                        aria-haspopup="true"
                        component={"h3"}>
                        <b>{dataset.name}</b>
                    </Link>
                    }
                    {dataset && <small style={{whiteSpace: 'nowrap'}}>&nbsp;
                        {hasSelection && shape != null && intFormat(selection.size) + ' / '}
                        {shape != null && intFormat(shape[0]) + ' cells'}
                    </small>}

                    {dataset != null && <Tabs
                        value={tab}
                        indicatorColor="primary"
                        textColor="primary"
                        onChange={this.handleTabChange}
                    >
                        <AntTab value="embedding" label="Embeddings"/>
                        <AntTab value="distribution" label="Distributions" disabled={distributionData.length === 0}/>
                        <AntTab value="composition" label="Composition" disabled={obsCat.length < 2}/>
                        {jobResults.length > 0 && <AntTab value="results" label="Results"/>}
                    </Tabs>}


                    <div style={{marginLeft: 'auto', whiteSpace: 'nowrap', overflow: 'hidden'}}>
                        {/*{serverInfo.brand && <Typography variant="h5"*/}
                        {/*                                 style={{*/}
                        {/*                                     display: 'inline-block',*/}
                        {/*                                     paddingRight: 6*/}
                        {/*                                 }}>{serverInfo.brand}</Typography>}*/}

                        {serverInfo.brand &&
                        <ReactMarkdown options={{
                            overrides: REACT_MD_OVERRIDES, wrapper: 'span', createElement: (type, props, children) => {
                                props.display = 'inline';
                                props.gutterBottom = false;
                                const elem = React.createElement(type, props, children);
                                return elem;
                            }
                        }}

                                       children={serverInfo.brand}/>}


                        {!loadingApp.loading && !isSignedOut && <DatasetSelector onChange={this.handleDataset}/>}
                        {showMoreMenu && <Tooltip title={'More'}>
                            <IconButton aria-label="Menu" aria-haspopup="true"
                                        onClick={this.handleMoreMenuOpen}>
                                <MoreVertIcon/>
                            </IconButton>
                        </Tooltip>}
                        {showMoreMenu && <Menu id="more-menu"
                                               anchorEl={this.state.moreMenuAnchorEl}
                                               anchorOrigin={{
                                                   vertical: 'top',
                                                   horizontal: 'right',
                                               }}

                                               transformOrigin={{
                                                   vertical: 'top',
                                                   horizontal: 'right',
                                               }} open={this.state.moreMenuOpen}
                                               onClose={this.handleMoreMenuClose}>
                            {showNewDataset && <MenuItem onClick={this.handleImportDataset}>
                                New Dataset
                            </MenuItem>}

                            {showEditDeleteDataset && <MenuItem onClick={this.handleSettings}>Edit Dataset</MenuItem>}
                            {showEditDeleteDataset && <MenuItem onClick={this.handleDelete}>Delete Dataset</MenuItem>}
                            {(showNewDataset || showEditDeleteDataset) && dataset != null && <Divider/>}
                            {dataset != null && <MenuItem onClick={this.copyLink}>Copy Link </MenuItem>}
                        </Menu>}

                        {<Tooltip title={"Toggle Light/Dark Theme"}>
                            <IconButton edge={false} className={this.props.chartOptions.darkMode ? 'cirro-active' : ''}
                                        aria-label="Toggle Theme" onClick={() => this.onDarkMode()}>
                                <Brightness2Icon/>
                            </IconButton>
                        </Tooltip>}
                        {dataset != null && <Tooltip title={'Help'}>
                            <IconButton aria-label="Help"
                                        onClick={this.handleHelp}>
                                <HelpIcon/>
                            </IconButton>
                        </Tooltip>}
                        {email != null && email !== '' &&
                        <Tooltip title={email}>
                            <IconButton aria-label="Menu" aria-haspopup="true"
                                        onClick={this.handleUserMenuOpen}>
                                <AccountCircle/>
                            </IconButton>
                        </Tooltip>}
                        {email != null && email !== '' &&
                        <Menu id="menu-user"
                              anchorEl={this.state.userMenuAnchorEl}
                              anchorOrigin={{
                                  vertical: 'top',
                                  horizontal: 'right',
                              }}

                              transformOrigin={{
                                  vertical: 'top',
                                  horizontal: 'right',
                              }} open={this.state.userMenuOpen}
                              onClose={this.handleUserMenuClose}>
                            <MenuItem onClick={this.handleLogout}>Sign Out</MenuItem>
                        </Menu>}
                        {isSignedOut && <Button style={{whiteSpace: 'nowrap'}} variant="outlined" color="primary"
                                                onClick={this.props.handleLogin}>Sign In</Button>}
                    </div>
                </Toolbar>
            </AppBar>

        );
    }
}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        binSummary: state.binSummary,
        binValues: state.binValues,
        chartOptions: state.chartOptions,
        combineDatasetFilters: state.combineDatasetFilters,
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
        datasetFilter: state.datasetFilter,
        dialog: state.dialog,
        distributionData: state.distributionData,
        distributionPlotInterpolator: state.distributionPlotInterpolator,
        email: state.email,
        embeddingLabels: state.embeddingLabels,
        embeddings: state.embeddings,
        interpolator: state.interpolator,
        jobResults: state.jobResults,
        loading: state.loading,
        loadingApp: state.loadingApp,
        markerOpacity: state.markerOpacity,
        message: state.message,
        pointSize: state.pointSize,
        savedDatasetState: state.savedDatasetState,
        searchTokens: state.searchTokens,
        selection: state.selection,
        serverInfo: state.serverInfo,
        tab: state.tab,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        user: state.user
    };
};
const mapDispatchToProps = (dispatch, ownProps) => {
    return {
        handleTab: (value) => {
            dispatch(setTab(value));
        },
        setMessage: (value) => {
            dispatch(setMessage(value));
        },
        handleLogin: () => {
            dispatch(login());
        },
        handleLogout: () => {
            dispatch(logout());
        },

        handleSavedDatasetState: value => {
            dispatch(setSavedDatasetState(value));
        },

        handleDataset: value => {
            dispatch(setDataset(value));
        },
        handleDialog: (value) => {
            dispatch(setDialog(value));
        },
        handleChartOptions: (value) => {
            dispatch(setChartOptions(value));
        }
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(AppHeader));



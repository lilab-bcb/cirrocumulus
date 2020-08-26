import {Divider, IconButton, Menu, Tooltip} from '@material-ui/core';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import {withStyles} from '@material-ui/core/styles';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import Toolbar from '@material-ui/core/Toolbar';
import AccountCircle from '@material-ui/icons/AccountCircle';
import HelpIcon from '@material-ui/icons/Help';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import React from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    login,
    logout,
    setDataset,
    setDialog,
    setMessage,
    setTab,
} from './actions';
import {drawerWidth} from './App';
import {intFormat} from './formatters';
import {
    DEFAULT_DARK_MODE,
    DEFAULT_INTERPOLATOR,
    DEFAULT_MARKER_OPACITY,
    DEFAULT_SHOW_AXIS,
    DEFAULT_SHOW_FOG,
    DEFAULT_SHOW_LABELS,
    DEFAULT_UNSELECTED_MARKER_OPACITY
} from "./reducers";


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

    copyLink = (event) => {

        const {chartOptions, dataset, embeddings, searchTokens, datasetFilter, interpolator, markerOpacity, unselectedMarkerOpacity, dotPlotData} = this.props;
        let linkText = window.location.protocol + '//' + window.location.host;

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
        let jsonChartOptions = {};

       const defaultChartOptions = {showLabels: DEFAULT_SHOW_LABELS, showAxis: DEFAULT_SHOW_AXIS,
           showFog: DEFAULT_SHOW_FOG, darkMode: DEFAULT_DARK_MODE};

        for (let key in defaultChartOptions) {
            let value = chartOptions[key];
            if (value !== defaultChartOptions[key]) {
                jsonChartOptions[key] = value;
            }
        }
        if (Object.keys(jsonChartOptions).length > 0) {
            json.chartOptions = jsonChartOptions;
        }

        if (searchTokens.length > 0) {
            json.q = searchTokens;
        }

        let datasetFilterJson = {};
        for (let key in datasetFilter) {
            let value = datasetFilter[key];
            if (window.Array.isArray(value)) {
                datasetFilterJson[key] = value;
            } else {
                if (value.operation !== '' && !isNaN(value.value) && value.value != null) {
                    datasetFilterJson[key] = {operation: value.operation, value: value.value};
                }
            }
        }
        if (Object.keys(datasetFilterJson).length > 0) {
            json.datasetFilter = datasetFilterJson;
        }

        if (markerOpacity !== DEFAULT_MARKER_OPACITY) {
            json.markerOpacity = markerOpacity;
        }
        if (unselectedMarkerOpacity !== DEFAULT_UNSELECTED_MARKER_OPACITY) {
            json.unselectedMarkerOpacity = unselectedMarkerOpacity;
        }

        if (dotPlotData && dotPlotData.length > 0) {
            let sortOrder = {};
            dotPlotData.forEach(data => {
                sortOrder[data.name] = data.sortBy;
            });
            json.sort = sortOrder;
        }
        if (interpolator.name !== DEFAULT_INTERPOLATOR) {
            json.colorScheme = interpolator.name;
        }
        linkText += '?q=' + JSON.stringify(json);
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


    handleLogout = () => {
        this.setState({userMenuOpen: false});
        this.props.handleLogout();
    };

    handleImportDataset = (event) => {
        this.props.handleDialog(IMPORT_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };
    handleDataset = (event) => {
        this.props.handleDataset(event.target.value);
    };
    handleSettings = (event) => {
        this.props.handleDialog(EDIT_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };

    handleDelete = (event) => {
        this.props.handleDialog(DELETE_DATASET_DIALOG);
        this.setState({moreMenuOpen: false});
    };

    render() {
        const {
            dataset, loadingApp, email, datasetChoices, selection, classes, serverInfo, tab, user
        } = this.props;
        const shape = dataset != null && dataset.shape != null ? dataset.shape : [0, 0];
        const hasSelection = dataset != null && shape[0] > 0 && !isNaN(selection.count);
        const showNumberOfCells = !hasSelection && dataset != null && !(selection.count > 0) && shape[0] > 0 && (selection.count !== shape[0]);


        const showNewDataset = user != null && user.importer;
        const showEditDeleteDataset = dataset !== null && dataset.owner;
        const showMoreMenu = showNewDataset || dataset != null;

        return (

            <AppBar position="fixed" color="default" className={classes.appBar}>
                <Toolbar variant="dense">


                    {datasetChoices.length > 0 && <Select
                        style={{marginRight: 6}}
                        disableUnderline={true}
                        displayEmpty={true}
                        value={dataset == null ? '' : dataset.id}
                        onChange={this.handleDataset}
                        inputProps={{
                            name: 'dataset',
                            id: 'dataset-id',
                        }}
                    >
                        <MenuItem key="" value="" disabled>
                            Choose a dataset
                        </MenuItem>
                        {datasetChoices.map(dataset => <MenuItem
                            key={dataset.id} value={dataset.id}>{dataset.name}</MenuItem>)}
                    </Select>}

                    <div className={"cirro-condensed"} style={{display: 'inline-block'}}>
                        {hasSelection && intFormat(selection.count)}
                        {hasSelection && ' / ' + intFormat(shape[0]) + ' cells'}
                        {showNumberOfCells && intFormat(shape[0]) + ' cells'}
                    </div>

                    {dataset != null && <Tabs
                        value={tab}
                        indicatorColor="primary"
                        textColor="primary"
                        onChange={this.handleTabChange}
                        aria-label="view"

                    >
                        <AntTab value="embedding" label="Embeddings"/>
                        <AntTab value="dot_plot" label="Dot Plot"/>
                    </Tabs>}

                    <div className={"cirro-condensed"}>
                        {/*<CloudIcon style={{verticalAlign: 'bottom'}} fontSize={'large'}/>*/}
                        {/*<h3*/}
                        {/*    style={{display: 'inline', marginRight: 20}}>Cirro</h3>*/}


                    </div>
                    <div style={{marginLeft: 'auto'}}>

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
                        <Tooltip title={'Help'}>
                            <IconButton aria-label="Help"
                                        onClick={this.handleHelp}>
                                <HelpIcon/>
                            </IconButton>
                        </Tooltip>

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


                        {!loadingApp.loading && email == null && serverInfo.clientId !== '' &&
                        <Button style={{whiteSpace: 'nowrap'}} variant="outlined" color="primary"
                                onClick={this.props.handleLogin}>Sign In</Button>}
                    </div>
                </Toolbar>
            </AppBar>

        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        chartOptions: state.chartOptions,
        embeddings: state.embeddings,
        searchTokens: state.searchTokens,
        binSummary: state.binSummary,
        binValues: state.binValues,
        combineDatasetFilters: state.combineDatasetFilters,
        datasetFilter: state.datasetFilter,
        markerOpacity: state.markerOpacity,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,

        datasetChoices: state.datasetChoices,
        dialog: state.dialog,
        dotPlotData: state.dotPlotData,
        email: state.email,

        interpolator: state.interpolator,
        loading: state.loading,
        loadingApp: state.loadingApp,

        message: state.message,

        selection: state.selection,
        serverInfo: state.serverInfo,
        user: state.user,
        tab: state.tab
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
        handleDataset: value => {
            dispatch(setDataset(value));
        },
        handleDialog: (value) => {
            dispatch(setDialog(value));
        }
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(AppHeader));



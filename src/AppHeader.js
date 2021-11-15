import {Divider, IconButton, Menu, Tooltip} from '@mui/material';
import AppBar from '@mui/material/AppBar';
import Button from '@mui/material/Button';
import MenuItem from '@mui/material/MenuItem';
import Popover from '@mui/material/Popover';
import Tab from '@mui/material/Tab';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';
import AccountCircle from '@mui/icons-material/AccountCircle';
import Brightness2Icon from '@mui/icons-material/Brightness3';
import HelpIcon from '@mui/icons-material/Help';
import MoreVertIcon from '@mui/icons-material/MoreVert';
import ReactMarkdown from 'markdown-to-jsx';
import React, {useState} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    getDatasetStateJson,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    login,
    logout,
    setChartOptions,
    setDataset,
    setDialog,
    setDrawerOpen,
    setMessage,
    setSavedDatasetState,
    setTab
} from './actions';
import CirroIcon from './CirroIcon';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import DatasetSelector from './DatasetSelector';
import {intFormat} from './formatters';
import {
    copyToClipboard,
    FEATURE_TYPE,
    REACT_MD_OVERRIDES,
    SERVER_CAPABILITY_ADD_DATASET,
    SERVER_CAPABILITY_DELETE_DATASET,
    SERVER_CAPABILITY_EDIT_DATASET
} from './util';
import Box from '@mui/material/Box';
import Tabs from '@mui/material/Tabs';
import MenuIcon from '@mui/icons-material/Menu';
import {DATASET_FIELDS} from './EditNewDatasetDialog';


function AppHeader(props) {

    const {
        chartOptions,
        dataset,
        drawerOpen,
        distributionData,
        handleChartOptions,
        handleLogin,
        handleLogout,
        handleDataset,
        handleDrawerOpen,
        handleDialog,
        handleTab,
        loadingApp,
        handleSavedDatasetState,
        jobResults,
        email,
        savedDatasetState,
        selection,
        searchTokens,
        serverInfo,
        tab,
        user
    } = props;

    const [userMenuOpen, setUserMenuOpen] = useState(false);
    const [userMenuAnchorEl, setUserMenuAnchorEl] = useState(null);
    const [moreMenuOpen, setMoreMenuOpen] = useState(false);
    const [moreMenuAnchorEl, setMoreMenuAnchorEl] = useState(null);
    const [datasetDetailsEl, setDatasetDetailsEl] = useState(null);

    function onTabChange(event, value) {
        handleTab(value);
    }

    function toggleDrawer() {
        handleDrawerOpen(!drawerOpen);
    }

    function onUserMenuClose() {
        setUserMenuOpen(false);
    }

    function onMoreMenuClose() {
        setMoreMenuOpen(false);
    }

    function onHelp() {
        handleDialog(HELP_DIALOG);
    }

    function onUserMenuOpen(event) {
        setUserMenuOpen(true);
        setUserMenuAnchorEl(event.currentTarget);
    }

    function onMoreMenuOpen(event) {
        setMoreMenuOpen(true);
        setMoreMenuAnchorEl(event.currentTarget);
    }

    function getLinkJson() {
        return getDatasetStateJson(props);
    }

    function onDataset(id) {
        if (dataset != null) {
            const link = getLinkJson();
            link.dataset = null;
            savedDatasetState[dataset.id] = link;
            handleSavedDatasetState(savedDatasetState);
        }
        handleTab('embedding'); // embedding won't render unless visible
        handleDataset(id);
    }

    function copyLink(event) {
        let linkText = window.location.protocol + '//' + window.location.host + window.location.pathname;
        linkText += '#q=' + encodeURIComponent(JSON.stringify(getLinkJson()));
        copyToClipboard(linkText);
        setMessage('Link copied');
        setMoreMenuOpen(false);
    }


    function onDarkMode() {
        chartOptions.darkMode = !chartOptions.darkMode;
        handleChartOptions(chartOptions);
    }

    function onLogout() {
        setUserMenuOpen(false);
        handleLogout();
    }

    function onImportDataset(event) {
        handleDialog(IMPORT_DATASET_DIALOG);
        setMoreMenuOpen(false);
    }

    function onSettings(event) {
        handleDialog(EDIT_DATASET_DIALOG);
        setMoreMenuOpen(false);
    }

    function onDelete(event) {
        handleDialog(DELETE_DATASET_DIALOG);
        setMoreMenuOpen(false);
    }

    function onShowDatasetDetails(event) {
        setDatasetDetailsEl(event.currentTarget);
    }

    function onCloseDatasetDetails(event) {
        setDatasetDetailsEl(null);
    }


    const datasetDetailsOpen = Boolean(datasetDetailsEl);
    const shape = dataset != null && dataset.shape != null ? dataset.shape : null;
    const hasSelection = dataset != null && shape != null && shape[0] > 0 && selection != null;
    const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT).map(item => item.id);
    const showAddDataset = user != null && user.importer && !loadingApp.loading && serverInfo.capabilities.has(SERVER_CAPABILITY_ADD_DATASET);
    const showEditDataset = dataset !== null && dataset.owner && !loadingApp.loading && serverInfo.capabilities.has(SERVER_CAPABILITY_EDIT_DATASET);
    const showDeleteDataset = dataset !== null && dataset.owner && !loadingApp.loading && serverInfo.capabilities.has(SERVER_CAPABILITY_DELETE_DATASET);
    const showMoreMenu = (showAddDataset || showEditDataset || showDeleteDataset || dataset != null) && !loadingApp.loading;
    const isSignedOut = !loadingApp.loading && email == null && serverInfo.clientId !== '';


    return (

        <Box sx={{display: 'flex'}}><AppBar position="fixed" color={"default"}
                                            sx={{
                                                zIndex: (theme) => theme.zIndex.drawer + 1
                                            }}>
            {dataset != null && datasetDetailsOpen && <Popover
                id={"dataset-details"}
                open={datasetDetailsOpen}
                anchorEl={datasetDetailsEl}
                onClose={onCloseDatasetDetails}
                anchorOrigin={{
                    vertical: 'bottom',
                    horizontal: 'center'
                }}
                transformOrigin={{
                    vertical: 'top',
                    horizontal: 'center'
                }}
            >
                <Box style={{width: 500, padding: '1em'}}>
                    <Typography variant="h6">{dataset.name}</Typography>
                    {DATASET_FIELDS.filter(item => dataset[item.fieldName]).map(item => <div
                        key={item.fieldName}><Divider/>
                        <Typography variant={"subtitle2"}>{item.label}</Typography>{item.fieldName !== 'description' &&
                        <Typography variant="body2"> {dataset[item.fieldName]}</Typography>}
                        {item.fieldName === 'description' &&
                        <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                       children={dataset[item.fieldName]}/>}</div>)}
                </Box>
            </Popover>
            }
            <Toolbar variant="dense" style={{paddingLeft: 0}}>
                <IconButton
                    size="large"
                    color="inherit"
                    aria-label="toggle drawer"
                    onClick={toggleDrawer}
                >
                    <MenuIcon/>
                </IconButton>
                <CirroIcon/>
                <Typography variant="h5" sx={{paddingRight: 1}}>Cirro</Typography>

                <Typography sx={{paddingRight: 1}}
                            variant="subtitle2">&nbsp;{hasSelection && shape != null && intFormat(selection.size) + ' / '}
                    {shape != null && intFormat(shape[0]) + ' cells'}</Typography>
                {dataset && <IconButton
                    aria-label="Info"
                    onClick={onShowDatasetDetails}
                    aria-owns={datasetDetailsOpen ? 'dataset-details' : undefined}
                    aria-haspopup="true"
                    size="small">
                    <InfoOutlinedIcon/>
                </IconButton>}
                <Typography variant="subtitle2">{dataset != null ? dataset.name : ''}</Typography>
                <div style={{display: 'flex', marginLeft: 'auto'}}>
                    <Tabs textColor="inherit" indicatorColor="secondary" value={tab} onChange={onTabChange}>
                        <Tab data-testid="embedding-tab" value="embedding" label="Embeddings"
                             disabled={dataset == null}/>
                        <Tab data-testid="distributions-tab" value="distribution" label="Distributions"
                             disabled={dataset == null || distributionData.length === 0}/>
                        <Tab data-testid="composition-tab" value="composition" label="Composition"
                             disabled={dataset == null || obsCat.length < 2}/>
                        {<Tab data-testid="results-tab" value="results" label="Results"
                              disabled={dataset == null || jobResults.length === 0}/>}
                    </Tabs>
                    {serverInfo.brand &&
                    <ReactMarkdown options={{
                        overrides: REACT_MD_OVERRIDES, wrapper: 'span', createElement: (type, props, children) => {
                            props.display = 'inline';
                            props.gutterBottom = false;
                            return React.createElement(type, props, children);
                        }
                    }} children={serverInfo.brand}/>}
                    {!loadingApp.loading && !isSignedOut && <DatasetSelector onChange={onDataset}/>}
                    {showMoreMenu && <Tooltip title={'More'}>
                        <IconButton
                            aria-label="Menu"
                            aria-haspopup="true"
                            onClick={onMoreMenuOpen}
                            size="large">
                            <MoreVertIcon/>
                        </IconButton>
                    </Tooltip>}
                    {showMoreMenu && <Menu id="more-menu"
                                           anchorEl={moreMenuAnchorEl}
                                           anchorOrigin={{
                                               vertical: 'top',
                                               horizontal: 'right'
                                           }}
                                           transformOrigin={{
                                               vertical: 'top',
                                               horizontal: 'right'
                                           }} open={moreMenuOpen}
                                           onClose={onMoreMenuClose}>
                        {showAddDataset && <MenuItem onClick={onImportDataset}>
                            New Dataset
                        </MenuItem>}

                        {showEditDataset && <MenuItem onClick={onSettings}>Edit Dataset</MenuItem>}
                        {showDeleteDataset && <MenuItem onClick={onDelete}>Delete Dataset</MenuItem>}
                        {(showAddDataset || showEditDataset || showDeleteDataset) && dataset != null && <Divider/>}
                        {dataset != null && <MenuItem onClick={copyLink}>Copy Link </MenuItem>}
                    </Menu>}
                    {<Tooltip title={"Toggle Light/Dark Theme"}>
                        <IconButton
                            edge={false}
                            className={chartOptions.darkMode ? 'cirro-active' : ''}
                            aria-label="Toggle Theme"
                            onClick={() => onDarkMode()}
                            size="large">
                            <Brightness2Icon/>
                        </IconButton>
                    </Tooltip>}
                    {dataset != null && <Tooltip title={'Help'}>
                        <IconButton aria-label="Help" onClick={onHelp} size="large">
                            <HelpIcon/>
                        </IconButton>
                    </Tooltip>}
                    {email != null && email !== '' &&
                    <Tooltip title={email}>
                        <IconButton
                            aria-label="Menu"
                            aria-haspopup="true"
                            onClick={onUserMenuOpen}
                            size="large">
                            <AccountCircle/>
                        </IconButton>
                    </Tooltip>}
                    {email != null && email !== '' &&
                    <Menu id="menu-user"
                          anchorEl={userMenuAnchorEl}
                          anchorOrigin={{
                              vertical: 'top',
                              horizontal: 'right'
                          }}

                          transformOrigin={{
                              vertical: 'top',
                              horizontal: 'right'
                          }} open={userMenuOpen}
                          onClose={onUserMenuClose}>
                        <MenuItem onClick={onLogout}>Sign Out</MenuItem>
                    </Menu>}
                    {isSignedOut && <Button style={{whiteSpace: 'nowrap'}} color="inherit"
                                            onClick={handleLogin}>Sign In</Button>}
                </div>
            </Toolbar>
        </AppBar>
        </Box>
    )
        ;

}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        chartOptions: state.chartOptions,
        combineDatasetFilters: state.combineDatasetFilters,
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
        datasetFilter: state.datasetFilter,
        datasetSelectorColumns: state.serverInfo.datasetSelectorColumns,
        dialog: state.dialog,
        distributionData: state.distributionData,
        distributionPlotInterpolator: state.distributionPlotInterpolator,
        drawerOpen: state.panel.drawerOpen,
        email: state.email,
        embeddingLabels: state.embeddingLabels,
        embeddings: state.embeddings,
        interpolator: state.interpolator,
        jobResults: state.jobResults,
        jobResultId: state.jobResultId,
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
const mapDispatchToProps = (dispatch) => {
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
        },
        handleDrawerOpen: (value) => {
            dispatch(setDrawerOpen(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(AppHeader));



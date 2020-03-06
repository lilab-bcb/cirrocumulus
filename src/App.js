import {IconButton, Menu, Snackbar, Tooltip} from '@material-ui/core';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import Chip from '@material-ui/core/Chip';
import CircularProgress from '@material-ui/core/CircularProgress';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import Drawer from '@material-ui/core/Drawer';
import Grid from '@material-ui/core/Grid';
import LinearProgress from '@material-ui/core/LinearProgress';
import Link from '@material-ui/core/Link';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import {withStyles} from '@material-ui/core/styles';
import Switch from '@material-ui/core/Switch';
import Toolbar from '@material-ui/core/Toolbar';
import AccountCircle from '@material-ui/icons/AccountCircle';
import CloseIcon from '@material-ui/icons/Close';
import CloudIcon from '@material-ui/icons/Cloud';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import React, {PureComponent} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    downloadSelectedIds,
    EDIT_DATASET_DIALOG,
    getDatasetFilterArray,
    IMPORT_DATASET_DIALOG,
    login,
    logout,
    removeDatasetFilter,
    SAVE_DATASET_FILTER_DIALOG,
    setCombineDatasetFilters,
    setDataset,
    setDialog,
    setMessage,
} from './actions';
import DeleteDatasetDialog from './DeleteDatasetDialog';
import DotPlotsPlotly from './DotPlots';
import EditDatasetDialog from './EditDatasetDialog';
import EmbeddingChartPlotly from './EmbeddingCharts';
import EmbedForm from './EmbedForm';
import {intFormat} from './formatters';
import {DEFAULT_INTERPOLATOR, DEFAULT_MARKER_OPACITY, DEFAULT_MARKER_SIZE} from "./reducers";
import SaveDatasetFilterDialog from './SaveDatasetFilterDialog';


const drawerWidth = 240;

const styles = theme => ({
    root: {
        display: 'flex',
    },
    appBar: {
        width: `calc(100% - ${drawerWidth}px)`,
        marginLeft: drawerWidth,
    },
    drawer: {
        width: drawerWidth,
        flexShrink: 0,
    },
    drawerPaper: {
        width: drawerWidth,
    },
    toolbar: theme.mixins.toolbar,
    content: {
        flexGrow: 1,
        backgroundColor: 'white',
        paddingTop: theme.spacing(6.5),
        paddingLeft: theme.spacing(1)
    },
});


class App extends PureComponent {


    constructor(props) {
        super(props);
        this.state = {
            userMenuOpen: false,
            userMenuAnchorEl: null,
            moreMenuOpen: false,
            moreMenuAnchorEl: null,

        };

    }

    handleUserMenuClose = () => {
        this.setState({userMenuOpen: false});
    };

    handleMoreMenuClose = () => {
        this.setState({moreMenuOpen: false});
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

    handleMessageClose = () => {
        this.props.setMessage(null);
    };
    handleUserMenuOpen = (event) => {
        this.setState({userMenuOpen: true, userMenuAnchorEl: event.currentTarget});
    };
    handleMoreMenuOpen = (event) => {
        this.setState({moreMenuOpen: true, moreMenuAnchorEl: event.currentTarget});
    };

    copyLink = () => {
        let linkText = window.location.protocol + '//' + window.location.host;

        let json = {
            dataset: this.props.dataset.id,
            embeddings: this.props.embeddings
        };
        if (this.props.features.length > 0) {
            json.features = this.props.features;
        }
        if (this.props.groupBy.length > 0) {
            json.groupBy = this.props.groupBy;
        }

        let datasetFilter = {};
        for (let key in this.props.datasetFilter) {
            let value = this.props.datasetFilter[key];
            if (window.Array.isArray(value)) {
                datasetFilter[key] = value;
            } else {
                if (value.operation !== '' && !isNaN(value.value) && value.value != null) {
                    datasetFilter[key] = {operation: value.operation, value: value.value};
                }
            }
        }
        json.datasetFilter = datasetFilter;
        if (this.props.markerSize !== DEFAULT_MARKER_SIZE) {
            json.markerSize = this.props.markerSize;
        }
        if (this.props.markerOpacity !== DEFAULT_MARKER_OPACITY) {
            json.markerOpacity = this.props.markerOpacity;
        }

        if (this.props.dotPlotData && this.props.dotPlotData.length > 0) {
            let sortOrder = {};
            this.props.dotPlotData.forEach(data => {
                sortOrder[data.name] = data.sortBy;
            });
            json.sort = sortOrder;
        }
        if (this.props.interpolator.name !== DEFAULT_INTERPOLATOR) {
            json.colorScheme = this.props.interpolator.name;
        }
        linkText += '?q=' + JSON.stringify(json);
        const el = document.createElement('textarea');
        el.value = linkText;
        el.setAttribute('readonly', '');
        el.style.position = 'absolute';
        el.style.left = '-9999px';
        document.body.appendChild(el);
        el.select();
        el.focus();
        document.execCommand('copy');
        document.body.removeChild(el);
        this.props.setMessage('Link copied');
        this.setState({moreMenuOpen: false});
    };

    handleSelectedCellsClick = (event) => {
        event.preventDefault();
        this.props.downloadSelectedIds();
    };


    handleLogout = () => {
        this.setState({userMenuOpen: false});
        this.props.handleLogout();
    };

    handleImportDataset = (event) => {
        this.props.handleDialog(IMPORT_DATASET_DIALOG);
    };
    handleDataset = (event) => {
        this.props.handleDataset(event.target.value);
    };
    handleSettings = (event) => {
        this.props.handleDialog(EDIT_DATASET_DIALOG);
    };

    handleDelete = (event) => {
        this.props.handleDialog(DELETE_DATASET_DIALOG);
    };

    render() {

        // tabs: 1. embedding, 2. grouped table with kde per feature, dotplot
        // need to add filter, selection
        const {classes, serverInfo, combineDatasetFilters} = this.props;
        const isEdit = this.props.savedDatasetFilter != null;
        let datasetFilters = getDatasetFilterArray(this.props.datasetFilter);
        let datasetFilterKeys = [];
        let isBrushing = false;
        datasetFilters.forEach(f => {
            if (typeof f[0] === 'object') {
                isBrushing = true;
            } else {
                datasetFilterKeys.push(f[0]);
            }
        });
        if (isBrushing) {
            datasetFilterKeys.push('selection');
        }

        const shape = this.props.dataset != null && this.props.dataset.shape != null ? this.props.dataset.shape : [0, 0];
        const hasSelection = this.props.dataset != null && shape[0] > 0 && !isNaN(this.props.selection.count);
        const showNumberOfCells = !hasSelection && this.props.dataset != null && !(this.props.selection.count > 0) && shape[0] > 0 && (this.props.selection.count !== shape[0]);
        return (
            <div className={classes.root}>
                {(this.props.dialog === EDIT_DATASET_DIALOG || this.props.dialog === IMPORT_DATASET_DIALOG) &&
                <EditDatasetDialog/>}
                {this.props.dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
                {this.props.dialog === SAVE_DATASET_FILTER_DIALOG && <SaveDatasetFilterDialog/>}
                <AppBar position="fixed" color="default" className={classes.appBar}>
                    <Toolbar variant="dense">
                        <div>
                            <CloudIcon style={{verticalAlign: 'bottom'}} fontSize={'large'}/>
                            <h3
                                style={{display: 'inline', marginRight: 20}}>Cirro</h3>

                            {this.props.email != null &&
                            <Select
                                disableUnderline={true}
                                displayEmpty={true}
                                value={this.props.dataset == null ? '' : this.props.dataset.id}
                                onChange={this.handleDataset}
                                inputProps={{
                                    name: 'dataset',
                                    id: 'dataset-id',
                                }}
                            > {this.props.datasetChoices.length > 0 && this.props.datasetChoices.length !== 1 &&
                            <MenuItem key="" value="" disabled>
                                Choose a dataset
                            </MenuItem>}
                                {this.props.datasetChoices.map(dataset => <MenuItem
                                    key={dataset.id} value={dataset.id}>{dataset.name}</MenuItem>)}


                            </Select>}
                            <div style={{display: 'inline-block', marginLeft: '10px'}}>
                                {hasSelection && (<Link title="Download selected ids" href="#"
                                                        onClick={this.handleSelectedCellsClick}>{intFormat(this.props.selection.count)}</Link>)}
                                {hasSelection && ' / ' + intFormat(shape[0]) + ' cells'}
                                {showNumberOfCells && intFormat(shape[0]) + ' cells'}
                            </div>

                            <div style={{display: 'inline-block', marginLeft: '10px'}}>

                                {datasetFilterKeys.map(key => {
                                    return <Chip
                                        size="small"
                                        onDelete={() => {
                                            this.onDatasetFilterChipDeleted(key);
                                        }}
                                        style={{marginRight: 2}}
                                        key={key}
                                        label={key}
                                        variant={'outlined'}
                                    />;
                                })}


                                {datasetFilters.length > 0 &&
                                <div style={{display: 'inline-block', marginLeft: '10px'}}><Link title="Clear" href="#"
                                                                                                 onClick={this.onDatasetFilterCleared}>Clear</Link>
                                </div>}
                                {datasetFilters.length > 0 &&
                                <div style={{display: 'inline-block', marginLeft: '10px'}}>
                                    <Grid component="label" container alignItems="center" spacing={0}>
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
                                </div>}
                                {datasetFilters.length > 0 && serverInfo.canWrite &&
                                <div style={{display: 'inline-block', marginLeft: '10px'}}>
                                    <Link style={{
                                        textOverflow: 'ellipsis',
                                        overflow: 'hidden',
                                        maxWidth: 140,
                                        display: 'inline-block',
                                        whiteSpace: 'nowrap',
                                        verticalAlign: 'text-bottom'
                                    }} title={isEdit ? 'Edit ' + this.props.savedDatasetFilter.name : 'Save'} href="#"
                                          onClick={this.onDatasetFilterSaved}>{isEdit ? 'Edit ' + this.props.savedDatasetFilter.name : 'Save'}</Link>
                                </div>}
                            </div>
                        </div>
                        <div style={{marginLeft: 'auto'}}>
                            <Tooltip title={'More'}>
                                <IconButton style={{marginLeft: 50}} aria-label="Menu" aria-haspopup="true"
                                            onClick={this.handleMoreMenuOpen}>
                                    <MoreVertIcon/>
                                </IconButton>
                            </Tooltip>
                            <Menu id="more-menu"
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
                                {this.props.user.importer &&
                                <MenuItem onClick={this.handleImportDataset}>
                                    Import Dataset
                                </MenuItem>}
                                {this.props.dataset !== null && this.props.dataset.owner &&
                                <MenuItem onClick={this.handleDelete}>Delete Dataset</MenuItem>}
                                {this.props.dataset !== null && this.props.dataset.owner &&
                                <MenuItem onClick={this.handleSettings}>Edit Dataset</MenuItem>}
                                {this.props.dataset != null &&
                                <MenuItem onClick={this.copyLink}>Copy Link
                                </MenuItem>}
                            </Menu>


                            {this.props.email != null &&
                            <Tooltip title={this.props.email}>
                                <IconButton style={{marginLeft: 50}} aria-label="Menu" aria-haspopup="true"
                                            onClick={this.handleUserMenuOpen}>
                                    <AccountCircle/>
                                </IconButton>
                            </Tooltip>}
                            {this.props.email != null &&
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


                            {!this.props.loadingApp.loading && this.props.email == null && serverInfo.clientId !== '' &&
                            <Button style={{whiteSpace: 'nowrap'}} variant="outlined" color="primary"
                                    onClick={this.props.handleLogin}>Sign In</Button>}
                        </div>
                    </Toolbar>
                </AppBar>
                <Drawer
                    className={classes.drawer}
                    variant="permanent"
                    classes={{
                        paper: classes.drawerPaper,
                    }}
                    anchor="left"
                >
                    {this.props.dataset != null && <EmbedForm key={this.props.dataset.id}/>}
                </Drawer>
                <div/>
                <main className={classes.content}>
                    {this.props.loadingApp.loading &&
                    <div><h2>Loading<LinearProgress style={{width: '90%'}} variant="determinate"
                                                    value={this.props.loadingApp.progress}/></h2>
                    </div>}
                    {this.props.dataset != null && <EmbeddingChartPlotly/>}
                    {this.props.dataset != null && <DotPlotsPlotly/>}
                </main>

                {this.props.loading && <Dialog aria-labelledby="loading-dialog-title" open={true}>
                    <DialogTitle id="loading-dialog-title"><CircularProgress size={20}/> Loading...</DialogTitle>
                </Dialog>}


                {this.props.message != null && <Snackbar
                    anchorOrigin={{
                        vertical: 'bottom',
                        horizontal: 'left',
                    }}
                    ContentProps={{
                        'aria-describedby': 'message-id',
                    }}
                    onClose={this.handleMessageClose}
                    open={true}
                    autoHideDuration={6000}
                    action={[
                        <IconButton
                            key="close"
                            aria-label="Close"
                            color="inherit"
                            onClick={this.handleMessageClose}
                        >
                            <CloseIcon/>
                        </IconButton>,
                    ]}
                    message={<span id="message-id">{this.props.message instanceof Error
                        ? this.props.message.message
                        : this.props.message}</span>}
                />}
            </div>
        );
    }
}

const mapStateToProps = state => {
    return {
        binSummary: state.binSummary,
        binValues: state.binValues,
        combineDatasetFilters: state.combineDatasetFilters,
        datasetFilter: state.datasetFilter,
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
        dialog: state.dialog,
        dotPlotData: state.dotPlotData,
        email: state.email,
        embeddings: state.embeddings,
        features: state.features,
        groupBy: state.groupBy,
        interpolator: state.interpolator,
        loading: state.loading,
        loadingApp: state.loadingApp,
        markerOpacity: state.markerOpacity,
        markerSize: state.markerSize,
        message: state.message,
        numberOfBins: state.numberOfBins,
        savedDatasetFilter: state.savedDatasetFilter,
        selection: state.selection,
        serverInfo: state.serverInfo,
        user: state.user,
        view3d: state.view3d

    };
};
const mapDispatchToProps = dispatch => {
    return {
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
        },
        handleCombineDatasetFilters: (value) => {
            dispatch(setCombineDatasetFilters(value));
        },
        setMessage: (value) => {
            dispatch(setMessage(value));
        },
        downloadSelectedIds: () => {
            dispatch(downloadSelectedIds());
        },
        removeDatasetFilter: (filter) => {
            dispatch(removeDatasetFilter(filter));
        }

    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(App));

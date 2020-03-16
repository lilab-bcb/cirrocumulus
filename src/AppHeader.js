import {IconButton, Menu, Tooltip} from '@material-ui/core';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import Chip from '@material-ui/core/Chip';
import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import {withStyles} from '@material-ui/core/styles';
import Switch from '@material-ui/core/Switch';
import Toolbar from '@material-ui/core/Toolbar';
import AccountCircle from '@material-ui/icons/AccountCircle';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import React from 'react';
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
import {intFormat} from './formatters';
import {DEFAULT_INTERPOLATOR, DEFAULT_MARKER_OPACITY} from "./reducers";

const drawerWidth = 240;


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
        if (Object.keys(datasetFilter).length > 0) {
            json.datasetFilter = datasetFilter;
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
        this.setState({moreMenuOpen: false});
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
        const {
            dataset, loadingApp, email, datasetChoices, selection, classes, serverInfo, combineDatasetFilters, savedDatasetFilter,
            datasetFilter, user
        } = this.props;
        const shape = dataset != null && dataset.shape != null ? dataset.shape : [0, 0];
        const hasSelection = dataset != null && shape[0] > 0 && !isNaN(selection.count);
        const showNumberOfCells = !hasSelection && dataset != null && !(selection.count > 0) && shape[0] > 0 && (selection.count !== shape[0]);

        const isEdit = savedDatasetFilter != null;
        let datasetFilters = getDatasetFilterArray(datasetFilter);
        const datasetFilterKeys = [];
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


        return (
            <AppBar position="fixed" color="default" className={classes.appBar}>
                <Toolbar variant="dense">
                    {email != null &&
                    <Select
                        style={{marginRight: 6}}
                        disableUnderline={true}
                        displayEmpty={true}
                        value={dataset == null ? '' : dataset.id}
                        onChange={this.handleDataset}
                        inputProps={{
                            name: 'dataset',
                            id: 'dataset-id',
                        }}
                    > {datasetChoices.length > 0 && datasetChoices.length !== 1 &&
                    <MenuItem key="" value="" disabled>
                        Choose a dataset
                    </MenuItem>}
                        {datasetChoices.map(dataset => <MenuItem
                            key={dataset.id} value={dataset.id}>{dataset.name}</MenuItem>)}
                    </Select>}

                    <div className={"cirro-condensed"}>
                        {/*<CloudIcon style={{verticalAlign: 'bottom'}} fontSize={'large'}/>*/}
                        {/*<h3*/}
                        {/*    style={{display: 'inline', marginRight: 20}}>Cirro</h3>*/}


                        <div style={{display: 'inline-block'}}>
                            {hasSelection && (<Link title="Download selected ids" href="#"
                                                    onClick={this.handleSelectedCellsClick}>{intFormat(selection.count)}</Link>)}
                            {hasSelection && ' / ' + intFormat(shape[0]) + ' cells'}
                            {showNumberOfCells && intFormat(shape[0]) + ' cells'}
                        </div>

                        <div style={{display: 'inline-block', marginLeft: '10px'}}>
                            {datasetFilters.length > 0 && <div style={{display: 'inline-block'}}>FILTER:</div>}
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
                            <div style={{display: 'inline-block', marginLeft: '10px'}}>
                                <Button size="small" color={'primary'}
                                        onClick={this.onDatasetFilterCleared}>Clear</Button>
                                <Button size="small" color={'primary'}
                                        onClick={this.onDatasetFilterSaved}>Save</Button>
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
                            {user.importer &&
                            <MenuItem onClick={this.handleImportDataset}>
                                Import Dataset
                            </MenuItem>}
                            {dataset !== null && dataset.owner &&
                            <MenuItem onClick={this.handleDelete}>Delete Dataset</MenuItem>}
                            {dataset !== null && dataset.owner &&
                            <MenuItem onClick={this.handleSettings}>Edit Dataset</MenuItem>}
                            {dataset != null &&
                            <MenuItem onClick={this.copyLink}>Copy Link
                            </MenuItem>}
                        </Menu>


                        {email != null &&
                        <Tooltip title={email}>
                            <IconButton style={{marginLeft: 50}} aria-label="Menu" aria-haspopup="true"
                                        onClick={this.handleUserMenuOpen}>
                                <AccountCircle/>
                            </IconButton>
                        </Tooltip>}
                        {email != null &&
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
        features: state.features,
        embeddings: state.embeddings,
        groupBy: state.groupBy,
        binSummary: state.binSummary,
        binValues: state.binValues,
        combineDatasetFilters: state.combineDatasetFilters,
        datasetFilter: state.datasetFilter,

        datasetChoices: state.datasetChoices,
        dialog: state.dialog,
        dotPlotData: state.dotPlotData,
        email: state.email,

        interpolator: state.interpolator,
        loading: state.loading,
        loadingApp: state.loadingApp,
        markerOpacity: state.markerOpacity,

        message: state.message,
        numberOfBins: state.numberOfBins,
        savedDatasetFilter: state.savedDatasetFilter,
        selection: state.selection,
        serverInfo: state.serverInfo,
        user: state.user,
        view3d: state.view3d
    };
};
const mapDispatchToProps = (dispatch, ownProps) => {
    return {
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
        },
        handleCombineDatasetFilters: (value) => {
            dispatch(setCombineDatasetFilters(value));
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
)(AppHeader));



import {IconButton, Menu, Snackbar, Tooltip} from '@material-ui/core';
import AppBar from '@material-ui/core/AppBar';
import Button from '@material-ui/core/Button';
import CircularProgress from '@material-ui/core/CircularProgress';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import Drawer from '@material-ui/core/Drawer';
import Input from '@material-ui/core/Input';
import Link from '@material-ui/core/Link';
import MenuItem from '@material-ui/core/MenuItem';
import Popover from '@material-ui/core/Popover';
import Select from '@material-ui/core/Select';
import withStyles from "@material-ui/core/styles/withStyles";
import Toolbar from '@material-ui/core/Toolbar';
import AccountCircle from '@material-ui/icons/AccountCircle';
import CloseIcon from '@material-ui/icons/Close';
import CloudIcon from '@material-ui/icons/Cloud';
import DeleteIcon from '@material-ui/icons/Delete';
import LinkIcon from '@material-ui/icons/Link';
import SettingsIcon from '@material-ui/icons/Settings';
import React, {PureComponent} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    downloadSelectedIds,
    EDIT_DATASET_DIALOG,
    IMPORT_DATASET_DIALOG,
    login,
    logout,
    setDataset,
    setDialog,
    setMessage,
} from './actions';
import DeleteDatasetDialog from './DeleteDatasetDialog';
import DotPlotsPlotly from './DotPlotsPlotly';
import EditDatasetDialog from './EditDatasetDialog';
import EmbeddingChartPlotly from './EmbeddingChartPlotly';
import EmbedForm from './EmbedForm';
import {intFormat} from './formatters';
import {
    DEFAULT_BIN_SUMMARY,
    DEFAULT_INTERPOLATOR,
    DEFAULT_MARKER_OPACITY,
    DEFAULT_MARKER_SIZE,
    DEFAULT_NUMBER_BINS
} from "./reducers";


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
        backgroundColor: theme.palette.background.default,
        paddingTop: theme.spacing(6.5),
        paddingLeft: theme.spacing(1)
    },
});


class App extends PureComponent {


    constructor(props) {
        super(props);
        this.state = {
            userMenuOpen: false,
            linkMenuOpen: false,
            userMenuAnchorEl: null,
            linkMenuAnchorEl: null,
            linkText: null,
        };
        this.linkRef = React.createRef();

    }

    handleUserMenuClose = () => {
        this.setState({userMenuOpen: false});
    };

    handleMessageClose = () => {
        this.props.setMessage(null);
    };
    handleUserMenuOpen = (event) => {
        this.setState({userMenuOpen: true, userMenuAnchorEl: event.currentTarget});
    };

    handleLinkMenuClose = (event) => {
        this.setState({linkMenuOpen: false, linkMenuAnchorEl: null});
    };

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (this.state.linkMenuOpen) {
            window.requestAnimationFrame(() => {
                this.linkRef.current.focus();
                this.linkRef.current.select();
            });
        }
    }

    copyLink = () => {
        this.linkRef.current.focus();
        this.linkRef.current.select();
        document.execCommand('copy');
        this.props.setMessage('Link copied');
    };

    handleSelectedCellsClick = (event) => {
        event.preventDefault();
        this.props.downloadSelectedIds();
    };

    handleLinkMenuOpen = (event) => {
        let linkText = window.location.protocol + '//' + window.location.host;


        let json = {
            dataset: this.props.dataset.id,
            embedding: this.props.embedding
        };
        if (this.props.features.length > 0) {
            json.features = this.props.features;
        }
        if (this.props.groupBy.length > 0) {
            json.groupBy = this.props.groupBy;
        }
        let categoricalFilter = this.props.categoricalFilter;
        if (Object.keys(categoricalFilter).length > 0) {
            json.categoricalFilter = this.props.categoricalFilter;
        }

        if (this.props.markerSize !== DEFAULT_MARKER_SIZE) {
            json.markerSize = this.props.markerSize;
        }
        if (this.props.markerOpacity !== DEFAULT_MARKER_OPACITY) {
            json.markerOpacity = this.props.markerOpacity;
        }
        if (this.props.binValues) {
            json.binValues = true;
            if (this.props.binSummary !== DEFAULT_BIN_SUMMARY) {
                json.binSummary = this.props.binSummary;
            }
            if (this.props.numberOfBins !== DEFAULT_NUMBER_BINS) {
                json.numberOfBins = this.props.numberOfBins;
            }
        }
        if (this.props.view3d) {
            json['3d'] = true;
        }

        if (this.props.interpolator.name !== DEFAULT_INTERPOLATOR) {
            json.colorScheme = this.props.interpolator.name;
        }

        linkText += '?q=' + JSON.stringify(json);
        this.setState({linkMenuOpen: true, linkMenuAnchorEl: event.currentTarget, linkText: linkText});

    };

    handleLogout = () => {
        this.setState({userMenuOpen: false});
        this.props.handleLogout();
    };

    handleDataset = (event) => {
        if (event.target.value === 'importDataset') {
            this.props.handleDialog(IMPORT_DATASET_DIALOG);
        } else {
            this.props.handleDataset(event.target.value);
        }
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
        const {classes} = this.props;
        const hasSelection = this.props.dataset != null && this.props.dataset.nObs > 0 && !isNaN(this.props.selection.count);
        const showNumberOfCells = !hasSelection && this.props.dataset != null && !(this.props.selection.count > 0) && this.props.dataset.nObs > 0 && (this.props.selection.count !== this.props.dataset.nObs);
        return (
            <div className={classes.root}>

                {(this.props.dialog === EDIT_DATASET_DIALOG || this.props.dialog === IMPORT_DATASET_DIALOG) &&
                <EditDatasetDialog/>}
                {this.props.dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
                <AppBar position="fixed" color="default" className={classes.appBar}>
                    <Toolbar variant="dense">
                        <div>
                            <CloudIcon style={{verticalAlign: 'bottom'}} fontSize={'large'}/>
                            <h3
                                style={{display: 'inline', marginRight: 20}}>Cirro</h3>
                            {this.props.loadingApp && <span>Loading...</span>}
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
                            > {this.props.datasetChoices.length > 0 && <MenuItem key="" value="" disabled>
                                Choose a dataset
                            </MenuItem>}
                                {this.props.datasetChoices.map(dataset => <MenuItem
                                    key={dataset.id} value={dataset.id}>{dataset.name}</MenuItem>)}

                                {this.props.user.importer && this.props.datasetChoices.length > 0 && <hr/>}
                                {this.props.user.importer && <MenuItem key="importDataset" value="importDataset">
                                    Import...
                                </MenuItem>}
                            </Select>}
                            <div style={{display: 'inline-block', marginLeft: '10px'}}>
                                {hasSelection && (<Link title="Download selected ids" href="#"
                                                        onClick={this.handleSelectedCellsClick}>{intFormat(this.props.selection.count)}</Link>)}
                                {hasSelection && ' / ' + intFormat(this.props.dataset.nObs) + ' cells'}
                                {showNumberOfCells && intFormat(this.props.dataset.nObs) + ' cells'}
                            </div>
                        </div>
                        <div style={{marginLeft: 'auto'}}>
                            {this.props.dataset != null &&
                            <Tooltip title="Link"><IconButton
                                aria-owns={this.state.linkMenuOpen ? 'link-popper' : undefined}
                                aria-haspopup="true" onClick={this.handleLinkMenuOpen} aria-label="Link">
                                <LinkIcon/>
                            </IconButton></Tooltip>}

                            {this.props.dataset !== null && this.props.dataset.owner &&
                            <Tooltip title="Edit"><IconButton onClick={this.handleSettings} aria-label="Edit">
                                <SettingsIcon/>
                            </IconButton></Tooltip>}

                            {this.props.dataset !== null && this.props.dataset.owner &&
                            <Tooltip title="Delete"><IconButton onClick={this.handleDelete} aria-label="Delete">
                                <DeleteIcon/>
                            </IconButton></Tooltip>}

                            {this.props.dataset != null && <Popover
                                id="link-popper"
                                open={this.state.linkMenuOpen}
                                anchorEl={this.state.linkMenuAnchorEl}
                                onClose={this.handleLinkMenuClose}

                                anchorOrigin={{
                                    vertical: 'bottom',
                                    horizontal: 'right',
                                }}
                                transformOrigin={{
                                    vertical: 'top',
                                    horizontal: 'right',
                                }}
                            >
                                <div style={{padding: 10}}>
                                    <h4>Copy link to share your current view</h4>
                                    <Button size="small" variant="outlined" onClick={this.copyLink}>
                                        Copy
                                    </Button> <Input autoFocus={true} inputRef={this.linkRef} readOnly={true}
                                                     value={this.state.linkText}></Input>
                                </div>
                            </Popover>
                            }
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


                            {!this.props.loadingApp && this.props.email == null &&
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
        categoricalFilter: state.categoricalFilter,
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
        dialog: state.dialog,
        email: state.email,
        embedding: state.embedding,
        features: state.features,
        groupBy: state.groupBy,
        interpolator: state.interpolator,
        loading: state.loading,
        loadingApp: state.loadingApp,
        markerOpacity: state.markerOpacity,
        markerSize: state.markerSize,
        message: state.message,
        numberOfBins: state.numberOfBins,
        selection: state.selection,
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
        setMessage: (value) => {
            dispatch(setMessage(value));
        },
        downloadSelectedIds: () => {
            dispatch(downloadSelectedIds());
        },

    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(App));

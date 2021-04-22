import {IconButton, Snackbar} from '@material-ui/core';
import CircularProgress from '@material-ui/core/CircularProgress';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import Drawer from '@material-ui/core/Drawer';
import LinearProgress from '@material-ui/core/LinearProgress';
import {createMuiTheme, withStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import CloseIcon from '@material-ui/icons/Close';
import {ThemeProvider} from '@material-ui/styles';
import React, {PureComponent} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    MORE_OPTIONS_DIALOG,
    SAVE_DATASET_FILTER_DIALOG,
    SAVE_FEATURE_SET_DIALOG,
    setDialog,
    setMessage,
} from './actions';
import AppHeader from './AppHeader';
import CompositionPlots from './CompositionPlots';
import DeleteDatasetDialog from './DeleteDatasetDialog';
import DistributionPlots from './DistributionPlots';
import DraggableDivider from './DraggableDivider';
import EditNewDatasetDialog from './EditNewDatasetDialog';
import EmbeddingChart from './EmbeddingChart';
import GalleryCharts from './GalleryCharts';
import HelpDialog from './HelpDialog';
import JobResultsPanel from './JobResultsPanel';
import LandingPage from './LandingPage';
import MoreOptionsDialog from './MoreOptionsDialog';
import SaveDatasetFilterDialog from './SaveDatasetFilterDialog';
import SaveSetDialog from './SaveSetDialog';
import SideBar from './SideBar';

const lightTheme = createMuiTheme(
    {
        "palette": {
            "type": "light"
        }
    }
);

const darkTheme = createMuiTheme(
    {
        "palette": {
            "type": "dark"
        }
    }
);
export const drawerWidth = 240;


const styles = (theme) => {
    return {
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
            paddingTop: theme.spacing(6.5),
            paddingLeft: theme.spacing(1)
        }
    };
};

class App extends PureComponent {

    constructor(props) {
        super(props);
        this.tooltipElementRef = React.createRef();
        this.galleryRef = React.createRef();
    }

    handleMessageClose = () => {
        this.props.setMessage(null);
    };


    onGallery = () => {
        window.scrollTo(0, this.galleryRef.current.offsetTop);
    };

    setTooltip = (text) => {
        text = text === '' || text == null ? '&nbsp;' : text;
        this.tooltipElementRef.current.innerHTML = text;
    };

    render() {

        // tabs: 1. embedding, 2. grouped table with kde per feature, dotplot
        // need to add filter, selection
        const {classes, chartOptions, dataset, dialog, loading, loadingApp, message, tab} = this.props;
        const theme = !chartOptions.darkMode ? lightTheme : darkTheme;
        const color = theme.palette.primary.main;
        const bgcolor = chartOptions.darkMode ? 'black' : 'white';
        const footerBackground = chartOptions.darkMode ? '#212121' : '#f5f5f5';
        return (<ThemeProvider theme={theme}>
                <div className={classes.root}>
                    {(dialog === EDIT_DATASET_DIALOG || dialog === IMPORT_DATASET_DIALOG) &&
                    <EditNewDatasetDialog/>}
                    {dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
                    {dialog === SAVE_DATASET_FILTER_DIALOG && <SaveDatasetFilterDialog/>}
                    {dialog === HELP_DIALOG && <HelpDialog/>}
                    {dialog === MORE_OPTIONS_DIALOG && <MoreOptionsDialog/>}
                    {dialog === SAVE_FEATURE_SET_DIALOG && <SaveSetDialog/>}
                    <AppHeader/>
                    <Drawer
                        className={classes.drawer}
                        variant="permanent"
                        classes={{
                            paper: classes.drawerPaper,
                        }}
                        anchor="left"
                    >
                        {dataset != null && <SideBar key={dataset.id}/>}
                    </Drawer>

                    <main style={{backgroundColor: bgcolor, color: color, paddingBottom: 24}}
                          className={classes.content}>
                        {loadingApp.loading &&
                        <div><h2>Loading<LinearProgress style={{width: '90%'}} variant="determinate"
                                                        value={loadingApp.progress}/></h2>
                        </div>}

                        {dataset == null && !loading && !loadingApp.loading && <div><LandingPage/></div>}
                        {dataset != null && <>
                            <div
                                role="tabpanel"
                                hidden={tab !== 'embedding'}
                            >
                                <EmbeddingChart onGallery={this.onGallery} setTooltip={this.setTooltip}/>
                                <DraggableDivider/>
                                <div ref={this.galleryRef}>
                                    <GalleryCharts/>
                                </div>
                            </div>
                            <div
                                role="tabpanel"
                                hidden={tab !== 'distribution'}
                            >
                                <DistributionPlots setTooltip={this.setTooltip}/>
                            </div>
                            <div
                                role="tabpanel"
                                hidden={tab !== 'composition'}
                            >
                                <CompositionPlots/>
                            </div>
                            <div
                                role="tabpanel"
                                hidden={tab !== 'results'}
                            >
                                <JobResultsPanel setTooltip={this.setTooltip}/>
                            </div>
                            <Typography className="cirro-condensed" color="textPrimary" ref={this.tooltipElementRef}
                                        style={{
                                            position: 'fixed',
                                            background: footerBackground,
                                            width: '100%',
                                            bottom: 0,
                                            top: 'auto',
                                            marginBottom: 0,
                                            whiteSpace: 'nowrap',
                                            textOverflow: 'ellipsis'
                                        }}>&nbsp;</Typography>
                        </>}
                    </main>

                    {loading && <Dialog aria-labelledby="loading-dialog-title" open={true}>
                        <DialogTitle id="loading-dialog-title"><CircularProgress size={20}/> Loading...</DialogTitle>
                    </Dialog>}


                    {message != null && <Snackbar
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
                        message={<span id="message-id">{message instanceof Error
                            ? message.message
                            : message}</span>}
                    />}
                </div>
            </ThemeProvider>
        );
    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        dataset: state.dataset,
        dialog: state.dialog,
        loading: state.loading,
        loadingApp: state.loadingApp,
        message: state.message,
        tab: state.tab
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleDialog: (value) => {
            dispatch(setDialog(value));
        },

        setMessage: (value) => {
            dispatch(setMessage(value));
        }
    };
};

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(App));

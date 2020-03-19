import {IconButton, Snackbar} from '@material-ui/core';
import CircularProgress from '@material-ui/core/CircularProgress';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import Drawer from '@material-ui/core/Drawer';
import LinearProgress from '@material-ui/core/LinearProgress';
import {withStyles} from '@material-ui/core/styles';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import CloseIcon from '@material-ui/icons/Close';
import React, {PureComponent} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    IMPORT_DATASET_DIALOG,
    SAVE_DATASET_FILTER_DIALOG,
    setDialog,
    setMessage,
} from './actions';
import AppHeader from './AppHeader';
import DeleteDatasetDialog from './DeleteDatasetDialog';
import DotPlotsPlotly from './DotPlots';
import EditDatasetDialog from './EditDatasetDialog';
import EmbeddingCharts from './EmbeddingCharts';
import EmbedForm from './EmbedForm';
import GalleryCharts from './GalleryCharts';
import SaveDatasetFilterDialog from './SaveDatasetFilterDialog';

const drawerWidth = 240;

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
        this.state = {tab: 'embedding'};
    }

    handleMessageClose = () => {
        this.props.setMessage(null);
    };

    handleTabChange = (event, value) => {
        this.setState({tab: value});
    };


    render() {

        // tabs: 1. embedding, 2. grouped table with kde per feature, dotplot
        // need to add filter, selection
        const {classes, dataset, dialog, loading, loadingApp, message} = this.props;

        return (
            <div className={classes.root}>
                {(dialog === EDIT_DATASET_DIALOG || dialog === IMPORT_DATASET_DIALOG) &&
                <EditDatasetDialog/>}
                {dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
                {dialog === SAVE_DATASET_FILTER_DIALOG && <SaveDatasetFilterDialog/>}

                <AppHeader/>
                <Drawer
                    className={classes.drawer}
                    variant="permanent"
                    classes={{
                        paper: classes.drawerPaper,
                    }}
                    anchor="left"
                >
                    {dataset != null && <EmbedForm key={dataset.id}/>}
                </Drawer>

                <main className={classes.content}>
                    {loadingApp.loading &&
                    <div><h2>Loading<LinearProgress style={{width: '90%'}} variant="determinate"
                                                    value={loadingApp.progress}/></h2>
                    </div>}


                    {dataset != null && <Tabs
                        value={this.state.tab}
                        indicatorColor="primary"
                        textColor="primary"
                        onChange={this.handleTabChange}
                        aria-label="view"

                    >
                        <AntTab value="embedding" label="Embedding"/>
                        <AntTab value="gallery" label="Gallery"/>
                        <AntTab value="dot_plot" label="Dot Plot"/>
                    </Tabs>}

                    <div
                        role="tabpanel"
                        hidden={this.state.tab !== 'embedding'}
                    >
                        {dataset != null && <EmbeddingCharts/>}
                    </div>
                    <div
                        role="tabpanel"
                        hidden={this.state.tab !== 'gallery'}
                    >
                        {dataset != null && <GalleryCharts/>}
                    </div>
                    <div
                        role="tabpanel"
                        hidden={this.state.tab !== 'dot_plot'}
                    >
                        {dataset != null && <DotPlotsPlotly/>}
                    </div>

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
        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        dialog: state.dialog,
        loading: state.loading,
        loadingApp: state.loadingApp,
        message: state.message,
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

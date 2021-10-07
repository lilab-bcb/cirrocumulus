import {IconButton, Snackbar} from '@mui/material';
import CircularProgress from '@mui/material/CircularProgress';
import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import Drawer from '@mui/material/Drawer';
import LinearProgress from '@mui/material/LinearProgress';

import Typography from '@mui/material/Typography';
import CloseIcon from '@mui/icons-material/Close';

import React, {PureComponent} from 'react';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    SAVE_DATASET_FILTER_DIALOG,
    SAVE_FEATURE_SET_DIALOG,
    setDialog,
    setMessage
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

import LandingPage from './LandingPage';
import SaveDatasetFilterDialog from './SaveDatasetViewDialog';
import SaveSetDialog from './SaveSetDialog';
import SideBar from './SideBar';
import {COMPARE_ACTIONS} from './job_config';
import Box from '@mui/material/Box';
import Toolbar from '@mui/material/Toolbar';
import {withTheme} from '@emotion/react';
import JobResultPanel from './JobResultPanel';


export const drawerWidth = 240;


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
        const {theme, dataset, dialog, loading, loadingApp, message, tab} = this.props;
        const color = theme.palette.primary.main;

        const footerBackground = theme.palette.background.paper;
        return (
            <Box sx={{display: 'flex'}}>
                {(dialog === EDIT_DATASET_DIALOG || dialog === IMPORT_DATASET_DIALOG) &&
                <EditNewDatasetDialog/>}
                {dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
                {dialog === SAVE_DATASET_FILTER_DIALOG && <SaveDatasetFilterDialog/>}
                {dialog === HELP_DIALOG && <HelpDialog/>}
                {dialog === SAVE_FEATURE_SET_DIALOG && <SaveSetDialog/>}
                <AppHeader/>
                <Drawer
                    variant="permanent"
                    anchor="left"
                    sx={{
                        width: drawerWidth,
                        flexShrink: 0,
                        '& .MuiDrawer-paper': {
                            width: drawerWidth,
                            boxSizing: 'border-box'
                        }
                    }}
                >
                    {dataset != null && <SideBar key={dataset.id} compareActions={COMPARE_ACTIONS}/>}
                </Drawer>


                <Box scomponent="main"
                     sx={{flexGrow: 1, paddingBottom: 24, color: color, backgroundColor: footerBackground}}>
                    <Toolbar/>
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
                            {<EmbeddingChart onGallery={this.onGallery} setTooltip={this.setTooltip}/>}
                            <DraggableDivider/>
                            <div ref={this.galleryRef}>
                                {<GalleryCharts/>}
                            </div>
                        </div>
                        <div
                            role="tabpanel"
                            hidden={tab !== 'distribution'}
                        >
                            {<DistributionPlots setTooltip={this.setTooltip}/>}
                        </div>
                        <div
                            role="tabpanel"
                            hidden={tab !== 'composition'}
                        >
                            {<CompositionPlots/>}
                        </div>
                        <div
                            role="tabpanel"
                            hidden={tab !== 'results'}
                        >
                            {<JobResultPanel setTooltip={this.setTooltip}/>}
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
                </Box>

                {loading && <Dialog aria-labelledby="loading-dialog-title" open={true}>
                    <DialogTitle id="loading-dialog-title"><CircularProgress
                        size={20}/> Loading...</DialogTitle>
                </Dialog>}


                {message != null && <Snackbar
                    anchorOrigin={{
                        vertical: 'bottom',
                        horizontal: 'left'
                    }}
                    ContentProps={{
                        'aria-describedby': 'message-id'
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
                            size="large">
                            <CloseIcon/>
                        </IconButton>
                    ]}
                    message={<span id="message-id">{message instanceof Error
                        ? message.message
                        : message}</span>}
                />}
            </Box>
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

export default withTheme(connect(
    mapStateToProps, mapDispatchToProps
)(App));

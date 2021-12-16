import {IconButton, Snackbar} from '@mui/material';
import CircularProgress from '@mui/material/CircularProgress';
import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import Drawer from '@mui/material/Drawer';
import LinearProgress from '@mui/material/LinearProgress';
import * as React from 'react';
import {useRef} from 'react';
import Box from '@mui/material/Box';
import Toolbar from '@mui/material/Toolbar';
import CloseIcon from '@mui/icons-material/Close';
import {connect} from 'react-redux';
import {
    DELETE_DATASET_DIALOG,
    EDIT_DATASET_DIALOG,
    HELP_DIALOG,
    IMPORT_DATASET_DIALOG,
    SAVE_DATASET_FILTER_DIALOG,
    SAVE_FEATURE_SET_DIALOG,
    setDialog,
    setDrawerOpen,
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
import {withTheme} from '@emotion/react';
import JobResultPanel from './JobResultPanel';


export const drawerWidth = 240;


function App(props) {
    const galleryRef = useRef();
    const {drawerOpen, theme, dataset, dialog, loading, loadingApp, message, setMessage, tab} = props;

    function handleMessageClose() {
        setMessage(null);
    }


    function onGallery() {
        window.scrollTo(0, galleryRef.current.offsetTop);
    }


    // tabs: 1. embedding, 2. grouped table with kde per feature, dotplot
    // need to add filter, selection

    const color = theme.palette.primary.main;
    const footerBackground = theme.palette.background.paper;
    return (
        <Box sx={{display: 'flex', backgroundColor: footerBackground}}>
            {(dialog === EDIT_DATASET_DIALOG || dialog === IMPORT_DATASET_DIALOG) &&
                <EditNewDatasetDialog/>}
            {dialog === DELETE_DATASET_DIALOG && <DeleteDatasetDialog/>}
            {dialog === SAVE_DATASET_FILTER_DIALOG && <SaveDatasetFilterDialog/>}
            {dialog === HELP_DIALOG && <HelpDialog/>}
            {dialog === SAVE_FEATURE_SET_DIALOG && <SaveSetDialog/>}
            <AppHeader/>
            <Drawer
                open={drawerOpen}
                variant="persistent"
                sx={{
                    width: drawerOpen ? drawerWidth : null,
                    flexShrink: 0,
                    '& .MuiDrawer-paper': {
                        width: drawerOpen ? drawerWidth : null,
                        boxSizing: 'border-box'
                    }
                }}
            >
                <Toolbar/>
                {dataset != null && <SideBar key={dataset.id} compareActions={COMPARE_ACTIONS}/>}
            </Drawer>
            <Box scomponent="main"
                 sx={{flexGrow: 1, marginLeft: 1, paddingBottom: 24, color: color, backgroundColor: footerBackground}}>
                <Toolbar/>
                {loadingApp.loading &&
                    <div><h2>Loading<LinearProgress style={{width: '90%'}} variant="determinate"
                                                    value={loadingApp.progress}/></h2>
                    </div>}

                {dataset == null && tab === 'embedding' && !loading && !loadingApp.loading &&
                    <div><LandingPage/></div>}
                {<>
                    {dataset != null && <div
                        role="tabpanel"
                        hidden={tab !== 'embedding'}
                    >
                        {<EmbeddingChart onGallery={onGallery}/>}
                        <DraggableDivider/>
                        <div ref={galleryRef}>
                            {<GalleryCharts/>}
                        </div>
                    </div>}
                    {dataset != null && <div
                        role="tabpanel"
                        hidden={tab !== 'distribution'}
                    >
                        {<DistributionPlots/>}
                    </div>}
                    {dataset != null && <div
                        role="tabpanel"
                        hidden={tab !== 'composition'}
                    >
                        {<CompositionPlots/>}
                    </div>}
                    {dataset != null && <div
                        role="tabpanel"
                        hidden={tab !== 'results'}
                    >
                        {<JobResultPanel/>}
                    </div>}
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
                onClose={handleMessageClose}
                open={true}
                autoHideDuration={6000}
                action={[
                    <IconButton
                        key="close"
                        aria-label="Close"
                        color="inherit"
                        onClick={handleMessageClose}
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

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        drawerOpen: state.panel.drawerOpen,
        dialog: state.dialog,
        loading: state.tasks.length > 0,
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
        handleDrawerOpen: (value) => {
            dispatch(setDrawerOpen(value));
        },
        setMessage: (value) => {
            dispatch(setMessage(value));
        }
    };
};

export default withTheme(connect(
    mapStateToProps, mapDispatchToProps
)(App));

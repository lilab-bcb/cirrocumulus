import {
    Dialog,
    DialogActions,
    DialogContent,
    DialogContentText,
    DialogTitle,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow
} from '@mui/material';
import Button from '@mui/material/Button';
import IconButton from '@mui/material/IconButton';
import DeleteIcon from '@mui/icons-material/Delete';
import React, {useState} from 'react';
import {deleteJobResult, setJobResultId} from './actions';
import withStyles from '@mui/styles/withStyles';
import {connect} from 'react-redux';
import {COMPARE_ACTIONS} from './job_config';
import Grid from '@mui/material/Grid';
import CancelIcon from '@mui/icons-material/Cancel';

function JobResultsSelector(props) {
    const [showDialog, setShowDialog] = useState(false);
    const [browseJob, setBrowseJob] = useState(null);
    const {email, handleJobResult, handleDeleteJob, jobResults, jobResultId} = props;

    function onBrowseJobs(event) {
        setShowDialog(true);
    }

    function onCloseJobs() {
        setShowDialog(false);
    }

    function onSelectJob(id) {
        handleJobResult(id);
        setShowDialog(false);
    }

    function onDeleteJob(event, job) {
        event.stopPropagation();
        setBrowseJob(job);
    }

    function onDeleteJobOK(id) {
        handleDeleteJob(id);
        setBrowseJob(null);
    }

    function onDeleteJobCancel() {
        setBrowseJob(null);
    }

    function isShowJobStatus() {
        for (let i = 0; i < jobResults.length; i++) {
            const isPrecomputed = ('' + jobResults[i].id).startsWith('cirro-');
            if (!isPrecomputed) {
                return true;
            }
        }
        return false;
    }

    const showJobStatus = isShowJobStatus();
    const jobTypeToName = {};
    const isShowingJob = jobResultId != null;
    COMPARE_ACTIONS.forEach(action => jobTypeToName[action.jobType] = action.title);

    function getTable() {
        return <Table size="small" stickyHeader={true}>
            <TableHead>
                <TableRow>
                    <TableCell>Name</TableCell>
                    <TableCell>Type</TableCell>
                    {showJobStatus && <TableCell>Status</TableCell>}
                </TableRow>
            </TableHead>
            <TableBody>
                {jobResults.map(jobResult => {
                    let text = jobResult.name;
                    if (jobResult.title) {
                        text += ' - ' + jobResult.title;
                    }
                    const isPrecomputed = ('' + jobResult.id).startsWith('cirro-');
                    const isComplete = isPrecomputed || (jobResult.status === 'complete' || jobResult.status === 'error');

                    const status = isPrecomputed ? 'complete' : jobResult.status;
                    const isJobOwner = email == jobResult.email || (email === null && jobResult.email === '');
                    const jobType = jobTypeToName[jobResult.type];
                    // const date = isPrecomputed ? '' : jobResult.submitted;
                    const showDelete = isJobOwner && !isPrecomputed && isComplete;
                    const showCancel = isJobOwner && !isPrecomputed && !isComplete;
                    return (
                        <TableRow key={jobResult.id}
                                  hover
                                  selected={jobResult.id === jobResultId}
                                  disabled={!isComplete}
                                  onClick={isComplete ? (event) => onSelectJob(jobResult.id) : null}
                                  role="checkbox"
                                  tabIndex={-1}>
                            <TableCell>{text}</TableCell>
                            <TableCell>{jobType}</TableCell>
                            {showJobStatus && <TableCell>{status}
                                {showDelete &&
                                <IconButton
                                    edge="end"
                                    aria-label="delete"
                                    onClick={(event) => onDeleteJob(event, jobResult)}
                                    size="small">
                                    <DeleteIcon/>
                                </IconButton>}
                                {showCancel &&
                                <IconButton
                                    edge="end"
                                    aria-label="cancel"
                                    onClick={(event) => onDeleteJob(event, jobResult)}
                                    size="small">
                                    <CancelIcon/>
                                </IconButton>}
                            </TableCell>}
                        </TableRow>
                    );
                })}
            </TableBody>
        </Table>;
    }

    return <>
        {isShowingJob && <Grid container alignItems="center">
            <Button size={"small"} onClick={onBrowseJobs} variant="outlined"
                    color="primary">Browse All Results</Button>
        </Grid>}
        <Dialog
            open={browseJob != null}
            onClose={onDeleteJobCancel}
            aria-labelledby="confirm-dialog-title"
            aria-describedby="conform-dialog-description"
        >
            <DialogTitle id="confirm-dialog-title">Confirm Delete</DialogTitle>
            <DialogContent>
                <DialogContentText id="confirm-dialog-description">
                    Are you sure you want to delete {browseJob ? browseJob.name : ''}?
                </DialogContentText>
            </DialogContent>
            <DialogActions>
                <Button onClick={onDeleteJobCancel}>
                    Cancel
                </Button>
                <Button variant="contained" onClick={e => onDeleteJobOK(browseJob.id)}
                        color="primary"
                        autoFocus>
                    OK
                </Button>
            </DialogActions>
        </Dialog>
        {isShowingJob &&
        <Dialog maxWidth={"xl"} fullWidth={true} onClose={onCloseJobs} aria-labelledby="job-results-title"
                open={showDialog && jobResults.length > 0}>
            <DialogTitle id="job-results-title" onClose={onCloseJobs}>
                Results
            </DialogTitle>
            <DialogContent>
                {getTable()}
            </DialogContent>
        </Dialog>}
        {!isShowingJob && getTable()}
    </>;
}

const mapStateToProps = state => {
        return {
            email: state.email,
            jobResultId: state.jobResultId,
            jobResults: state.jobResults
        };
    }
;
const mapDispatchToProps = (dispatch) => {
        return {
            handleDeleteJob: (payload) => {
                dispatch(deleteJobResult(payload));
            },
            handleJobResult: (payload) => {
                dispatch(setJobResultId(payload));
            }
        };
    }
;


export default withStyles({})(connect(
    mapStateToProps, mapDispatchToProps
)(JobResultsSelector));

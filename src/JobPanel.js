import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import DialogContent from '@mui/material/DialogContent';
import TextField from '@mui/material/TextField';
import {DialogActions, Divider, Switch, Typography} from '@mui/material';
import Button from '@mui/material/Button';
import {SERVER_CAPABILITY_JOBS} from './util';
import Tooltip from '@mui/material/Tooltip';
import Grid from '@mui/material/Grid';
import ButtonGroup from '@mui/material/ButtonGroup';
import {intFormat} from './formatters';
import FormControl from '@mui/material/FormControl';
import AutocompleteVirtualized from './AutocompleteVirtualized';
import FontDownloadRoundedIcon from '@mui/icons-material/FontDownloadRounded';
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';
import React, {useEffect, useState} from 'react';
import {datasetFilterToJson, setDialog, setTab, submitJob} from './actions';
import {connect} from 'react-redux';
import {getAnnotationOptions} from './ExplorePanel';

function JobPanel(props) {
    const {
        compareActions,
        combineDatasetFilters,
        dataset,
        datasetFilter,
        handleSubmitJob,
        selection,
        serverInfo,
        tab
    } = props;

    const [compareGroups, setCompareGroups] = useState('selected');
    const [group1, setGroup1] = useState(null);
    const [group1Count, setGroup1Count] = useState(null);

    const [group2, setGroup2] = useState(null);
    const [group2Count, setGroup2Count] = useState(null);

    const [jobName, setJobName] = useState('');
    const [jobParams, setJobParams] = useState(null);
    const [compareMenu, setCompareMenu] = useState(null);

    const [compareCategories, setCompareCategories] = useState([]);

    const obsCatOptions = getAnnotationOptions([], dataset.obsCat);

    useEffect(() => {
        setGroup1(null);
        setGroup1Count(null);
        setGroup2(null);
        setGroup2Count(null);
        setCompareCategories([]);
        setJobParams(null);
    }, [dataset]);

    function onCompareGroups(event) {
        setCompareGroups(event.target.checked ? 'selected' : 'all');
    }

    function onJobNameChange(event) {
        setJobName(event.target.value);
    }

    function onSubmitJobCancel() {
        setJobName('');
        setJobParams(null);
    }

    function onSubmitJobOK() {
        jobParams.name = jobName;
        handleSubmitJob(jobParams);
        setJobName('');
        setJobParams(null);
    }

    function onSetGroup(groupNumber) {
        const g = datasetFilterToJson(dataset, datasetFilter, combineDatasetFilters);
        const count = selection.size;
        if (groupNumber === 1) {
            setGroup1(g);
            setGroup1Count(count);
        } else {
            setGroup2(g);
            setGroup2Count(count);
        }
    }

    function onCompareCategories(event, value) {
        let values = [];
        value.forEach(val => {
            if (val.text !== undefined) {
                values.push(val.text);
            } else {
                values.push(val);
            }
        });
        setCompareCategories(values);

    }

    function onSubmitJob(jobType, version) {
        setJobName('');
        const p = {
            type: jobType,
            params: {version: version}
        };
        if (compareGroups === 'selected') {
            p.params.filter = group1;
            p.params.filter2 = group2;
        } else {
            p.params.obs = compareCategories;
        }
        setJobParams(p);
    }


    return <>
        <Dialog
            open={jobParams != null}
            onClose={onSubmitJobCancel}
            aria-labelledby="submit-job-dialog-title"
            aria-describedby="submit-job-dialog-description"
        >
            <DialogTitle id="submit-job-dialog-title">Submit</DialogTitle>
            <DialogContent>
                {/*<DialogContentText id="submit-job-dialog-description">*/}
                {/*    Job Details*/}
                {/*</DialogContentText>*/}
                <TextField
                    size={"small"}
                    onChange={onJobNameChange}
                    value={jobName}
                    autoFocus
                    margin="dense"
                    label="Job Name"
                    type="text"
                    fullWidth
                />
            </DialogContent>
            <DialogActions>
                <Button onClick={onSubmitJobCancel}>
                    Cancel
                </Button>
                <Button disabled={jobName === ''} variant="contained" onClick={e => onSubmitJobOK(jobName)}
                        color="primary">
                    Submit
                </Button>
            </DialogActions>
        </Dialog>

        {serverInfo.capabilities.has(SERVER_CAPABILITY_JOBS) &&
        <div style={tab === 'embedding' ? null : {display: 'none'}}>
            <Divider/>
            <Typography gutterBottom={false} component={"h1"}
                        style={{textTransform: 'uppercase'}}>Compare</Typography>

            <Grid alignContent={"flex-start"} container alignItems="center"
                  spacing={0}>
                {/*<Grid item><InputLabel shrink={true}>Combine</InputLabel></Grid>*/}
                <Tooltip
                    title={"Compare all pairs in a category"}><Grid item>ALL PAIRS</Grid></Tooltip>
                <Grid item>
                    <Switch
                        size="small"
                        checked={compareGroups === 'selected'}
                        onChange={onCompareGroups}
                    />
                </Grid>
                <Tooltip
                    title={"Choose a specific pair"}><Grid item>SELECTED</Grid></Tooltip>
            </Grid>
            {compareGroups === 'selected' &&
            <ButtonGroup variant="outlined" disabled={selection.size === 0}>
                <Tooltip
                    title={'Group one' + (group1Count ? ' (' + intFormat(group1Count) + ' cells)' : '')}>
                    <Button size={"small"}
                            onClick={event => onSetGroup(1)}>1</Button>
                </Tooltip>
                <Tooltip
                    title={'Group two' + (group2Count ? ' (' + intFormat(group2Count) + ' cells)' : '')}>
                    <Button size={"small"}
                            onClick={event => onSetGroup(2)}>2</Button>
                </Tooltip>
            </ButtonGroup>}
            {compareGroups === 'all' && <FormControl sx={{display: 'block'}}>
                <AutocompleteVirtualized label={"Metadata Categories"}
                                         testId={'compare-meta-input'}
                                         options={obsCatOptions}
                                         value={compareCategories}
                                         getOptionLabel={(option) => option.text}
                                         getChipIcon={(option) => {
                                             return <FontDownloadRoundedIcon
                                                 style={{
                                                     marginLeft: 4,
                                                     marginTop: 0,
                                                     marginRight: 0,
                                                     marginBottom: 0
                                                 }}
                                                 className={"MuiChip-deleteIcon MuiChip-deleteIconSmall"}/>;
                                         }}
                                         getOptionSelected={(option, value) => option.id === value}
                                         onChange={onCompareCategories}/>
            </FormControl>}

            <Button size={"small"} variant={"outlined"}
                    endIcon={<ArrowDropDownIcon/>}
                    disabled={(compareGroups === 'selected' && (group1 == null || group2 == null)) || (compareGroups === 'all' && compareCategories.length === 0)}
                    onClick={event => setCompareMenu(event.currentTarget)}>
                GO</Button>
            <Menu
                variant={"menu"}
                id="compare-menu"
                anchorEl={compareMenu}
                open={Boolean(compareMenu)}
                onClose={event => setCompareMenu(null)}
            >
                {compareActions.map(action => <MenuItem title={action.title}
                                                        key={action.jobType}
                                                        onClick={event => {
                                                            onSubmitJob(action.jobType, action.version);
                                                            setCompareMenu(null);
                                                        }}>{action.title}</MenuItem>)}

            </Menu>
        </div>}
    </>;
}


const mapStateToProps = state => {
        return {
            combineDatasetFilters: state.combineDatasetFilters,
            dataset: state.dataset,
            datasetFilter: state.datasetFilter,
            selection: state.selection,
            serverInfo: state.serverInfo,
            tab: state.tab
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
        return {
            handleDialog: (value) => {
                dispatch(setDialog(value));
            },
            handleTab: (value) => {
                dispatch(setTab(value));
            },
            handleSubmitJob: value => {
                dispatch(submitJob(value));
            }
        };
    }
;

export default (connect(
    mapStateToProps, mapDispatchToProps
)(JobPanel));
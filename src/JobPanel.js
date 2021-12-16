import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import DialogContent from '@mui/material/DialogContent';
import TextField from '@mui/material/TextField';
import {DialogActions, Divider, InputLabel, Typography} from '@mui/material';
import Button from '@mui/material/Button';
import {SERVER_CAPABILITY_JOBS} from './util';
import Tooltip from '@mui/material/Tooltip';
import ButtonGroup from '@mui/material/ButtonGroup';
import {intFormat} from './formatters';
import FormControl from '@mui/material/FormControl';
import AutocompleteVirtualized from './AutocompleteVirtualized';
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';
import React, {useEffect, useState} from 'react';
import {datasetFilterToJson, setDialog, setTab, submitJob} from './actions';
import {connect} from 'react-redux';
import {getAnnotationOptions} from './ExplorePanel';
import Select from '@mui/material/Select';
import withStyles from '@mui/styles/withStyles';

const styles = theme => ({
    formControl: {
        display: 'block', minWidth: 216, maxWidth: 216, marginBottom: theme.spacing(1)
    }, select: {
        width: 216
    }
});

function JobPanel(props) {
    const {
        classes,
        compareActions,
        combineDatasetFilters,
        dataset,
        datasetFilter,
        handleSubmitJob,
        selection,
        serverInfo,
        tab
    } = props;

    const [comparisonsOption, setComparisonsOption] = useState('one_vs_rest');
    const [group1, setGroup1] = useState(null);
    const [group1Count, setGroup1Count] = useState(null);

    const [group2, setGroup2] = useState(null);
    const [group2Count, setGroup2Count] = useState(null);

    const [jobName, setJobName] = useState('');
    const [jobParams, setJobParams] = useState(null);
    const [compareMenu, setCompareMenu] = useState(null);

    const [compareCategories, setCompareCategories] = useState([]);

    const obsCatOptions = getAnnotationOptions([], dataset.obsCat);

    const comparisonOptions = [{label: 'All Pairs', value: 'pairs', title: 'Compare all pairs of groups'}, {
        label: 'One vs. Rest', value: 'one_vs_rest', title: 'Compare group to all other groups'
    }, {label: 'Selection', value: 'selection', title: 'Compare selected cells'}];

    useEffect(() => {
        setGroup1(null);
        setGroup1Count(null);
        setGroup2(null);
        setGroup2Count(null);
        setCompareCategories([]);
        setJobParams(null);
    }, [dataset]);

    function onCompareGroups(event) {
        setComparisonsOption(event.target.value);
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
            type: jobType, params: {version: version}
        };
        if (comparisonsOption === 'selection') {
            p.params.filter = group1;
            p.params.filter2 = group2;
            p.params.pairs = '1';
        } else {
            if (comparisonsOption === 'pairs') {
                p.params.pairs = '1';
            }
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

                <FormControl className={classes.formControl}>
                    <InputLabel htmlFor="comparison_select">Comparisons</InputLabel>
                    <Select
                        label={"Comparisons"}
                        labelId={"comparison_select"}
                        size={"small"}
                        className={classes.select}
                        onChange={onCompareGroups}
                        value={comparisonsOption}
                        multiple={false}>
                        {comparisonOptions.map(item => (
                            <MenuItem key={item.value} value={item.value}>{item.label}</MenuItem>))}
                    </Select>
                </FormControl>


                {comparisonsOption === 'selection' &&
                    <ButtonGroup variant="outlined" disabled={selection == null || selection.size === 0}>
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
                {comparisonsOption !== 'selection' && <FormControl sx={{display: 'block'}}>
                    <AutocompleteVirtualized label={"Metadata Categories"}
                                             testId={'compare-meta-input'}
                                             options={obsCatOptions}
                                             value={compareCategories}
                                             getOptionLabel={(option) => option.text}
                                             getOptionSelected={(option, value) => option.id === value}
                                             helperText={'Select two or more categories to combine categories using "AND"'}
                                             onChange={onCompareCategories}/>
                </FormControl>}

                <Button size={"small"} variant={"outlined"}
                        endIcon={<ArrowDropDownIcon/>}
                        disabled={(comparisonsOption === 'selection' && (group1 == null || group2 == null)) || (comparisonsOption !== 'selection' && compareCategories.length === 0)}
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
};
const mapDispatchToProps = (dispatch) => {
    return {
        handleDialog: (value) => {
            dispatch(setDialog(value));
        }, handleTab: (value) => {
            dispatch(setTab(value));
        }, handleSubmitJob: value => {
            dispatch(submitJob(value));
        }
    };
};

export default withStyles(styles)(connect(mapStateToProps, mapDispatchToProps)(JobPanel));
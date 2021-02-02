import {InputLabel, MenuItem, Select, Tooltip} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import Input from '@material-ui/core/Input';
import Slider from '@material-ui/core/Slider';
import withStyles from '@material-ui/core/styles/withStyles';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {debounce, find} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setJobResults} from './actions';
import {EditableColorScheme} from './EditableColorScheme';
import {EditableSizeLegend} from './EditableSizeLegend';
import {sortAndFilterJobResult, sortByGroup, updateJob, updateTopNJobResult} from './JobResultsPanel';
import {createColorScale} from './util';


const styles = theme => ({
    formControl: {
        // display: 'block',
        margin: theme.spacing(0, 1)
    }
});

class JobResultOptions extends React.PureComponent {

    constructor(props) {
        super(props);
        this.handleJobResultDebounced = debounce(this.handleJobResultDebounced, 500);
        this.sortAndFilterDebounced = debounce(this.sortAndFilterDebounced, 500);
        this.updateTopNJobResultDebounced = debounce(this.updateTopNJobResultDebounced, 500);
        this.state = {
            forceUpdate: false,
        };
    }

    getJobResult = (id = null) => {
        const {jobResults, jobResultId} = this.props;
        if (id == null) {
            id = jobResultId;
        }
        return find(jobResults, item => item.id === id);
    };

    updateTopNJobResultDebounced = () => {
        updateTopNJobResult(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    sortAndFilterDebounced = () => {
        sortAndFilterJobResult(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    handleJobResultDebounced = (id) => {
        updateJob(this.getJobResult(id));
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onSizeReversedChange = (value) => {
        const jobResult = this.getJobResult();
        jobResult.sizeScaleReversed = value;
        jobResult.sizeScale = null;
        updateJob(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onSortChange = (event) => {
        const jobResult = this.getJobResult();
        jobResult.sortByGroup = event.target.value;
        sortByGroup(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onOptions = (options) => {
        const jobResult = this.getJobResult();
        jobResult.options = Object.assign(jobResult.options, options);
        jobResult.colorScale = null;
        jobResult.sizeScale = null;
        this.handleJobResultDebounced(this.props.jobResultId);
    };

    onInterpolator = (value) => {
        const jobResult = this.getJobResult();
        jobResult.interpolator = value;
        jobResult.colorScale = createColorScale(jobResult.interpolator).domain(jobResult.colorScale.domain()).clamp(true);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onColorChanged = (event) => {
        const jobResult = this.getJobResult();
        jobResult.colorScale = null;
        jobResult.color = event.target.value;
        updateJob(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onSizeChanged = (event) => {
        const jobResult = this.getJobResult();
        jobResult.sizeScale = null;
        jobResult.size = event.target.value;
        updateJob(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onNtopChange = (event, value) => {
        const jobResult = this.getJobResult();
        jobResult.ntopUI = value;
        jobResult.ntop = value;
        this.setState((prevState, props) => ({
            forceUpdate: !prevState.forceUpdate
        }));
        this.updateTopNJobResultDebounced();
    };

    onByChanged = (event) => {
        const jobResult = this.getJobResult();
        jobResult.by = event.target.value;
        sortAndFilterJobResult(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onOperationChanged = (filter, value) => {
        filter[1] = value;
        sortAndFilterJobResult(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onValueChange = (filter, value) => {
        filter[2] = parseFloat(value);
        filter[3] = value;
        this.setState((prevState, props) => ({
            forceUpdate: !prevState.forceUpdate
        }));
        this.sortAndFilterDebounced();
    };

    render() {
        const {
            jobResultId,
            classes,
            textColor
        } = this.props;
        if (jobResultId == null) {
            return null;
        }
        const jobResult = this.getJobResult();
        if (jobResult == null) {
            return null;
        }
        const {
            rowFilters,
            color,
            size,
            fields,
            by,
            colorScale,
            ntopUI,
            interpolator,
            sizeScale
        } = jobResult;

        return <React.Fragment>

            <div style={{marginTop: 16}}>
                <Typography
                    component={"h2"}>Rank Features</Typography>
                <FormControl style={{display: 'block'}} className={classes.formControl}>
                    <InputLabel shrink={true}>By</InputLabel>
                    <Select
                        input={<Input size={"small"}/>}
                        onChange={this.onByChanged}
                        value={by}
                    >
                        {fields.map(item => (
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl>
                <InputLabel style={{marginLeft: 8, marginTop: 8}} shrink={true}>Number of Features</InputLabel>
                <Slider
                    min={5}
                    max={50}
                    step={5}
                    style={{marginLeft: 10, width: '86%'}}
                    valueLabelDisplay="auto"
                    value={ntopUI}
                    onChange={this.onNtopChange}/>
                <FormControl style={{display: 'block'}} className={classes.formControl}>
                    <InputLabel shrink={true}>Sort</InputLabel>
                    <Select
                        input={<Input size={"small"}/>}
                        onChange={this.onSortChange}
                        value={jobResult.sortByGroup}
                    >
                        {jobResult.groups.map(item => (
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl>
            </div>
            <div style={{marginTop: 8}}></div>
            <Tooltip title="Filters are applied separately per cluster"><Typography
                component={"h2"}>Filters</Typography></Tooltip>

            {rowFilters.map(filter => {
                // [field, op, val, uiValue]
                const id = 'job_result' + filter[0];
                return <div key={filter[0]} style={{paddingTop: 8}}>
                    <FormControl className={classes.formControl}>
                        <InputLabel>{filter[0]}</InputLabel>
                        <Select
                            labelId={id + '_label'}
                            id={id}
                            style={{marginRight: 6}}
                            value={filter[1]}
                            onChange={event => this.onOperationChanged(filter, event.target.value)}
                        >
                            <MenuItem value={""}></MenuItem>
                            <MenuItem value={">"}>{">"}</MenuItem>
                            <MenuItem value={"<"}>{"<"}</MenuItem>
                            <MenuItem value={"="}>{"="}</MenuItem>
                            <MenuItem value={">="}>{">="}</MenuItem>
                            <MenuItem value={"<="}>{"<="}</MenuItem>
                            <MenuItem value={"!="}>{"!="}</MenuItem>
                        </Select>
                    </FormControl>
                    <TextField
                        onChange={event => this.onValueChange(filter, event.target.value)} value={filter[3]}
                        style={{maxWidth: 60, verticalAlign: 'bottom'}}/>
                </div>;
            })}


            <Typography style={{marginTop: 16}}
                        component={"h2"}>View Options</Typography>

            <FormControl className={classes.formControl}>
                <InputLabel shrink={true}>Color</InputLabel>
                <Select
                    input={<Input size={"small"}/>}
                    onChange={this.onColorChanged}
                    value={color}
                >
                    {fields.map(item => (
                        <MenuItem key={item} value={item}>{item}</MenuItem>
                    ))}
                </Select>
            </FormControl>

            <EditableColorScheme colorScale={colorScale}
                                 textColor={textColor}
                                 interpolator={interpolator}
                                 onOptions={this.onOptions}
                                 onInterpolator={this.onInterpolator}/>
            <FormControl style={{marginTop: 16}} className={classes.formControl}>
                <InputLabel shrink={true}>Size</InputLabel>
                <Select
                    input={<Input size={"small"}/>}
                    onChange={this.onSizeChanged}
                    value={size}
                >
                    {fields.map(item => (
                        <MenuItem key={item} value={item}>{item}</MenuItem>
                    ))}
                </Select>
            </FormControl>
            <EditableSizeLegend sizeScale={sizeScale} textColor={textColor}
                                onOptions={this.onOptions} showReversed={true}
                                onReversedChange={this.onSizeReversedChange}/>
        </React.Fragment>;
    }
}

const mapStateToProps = state => {
        return {
            textColor: state.chartOptions.darkMode ? 'white' : 'black',
            jobResults: state.jobResults,
            jobResultId: state.jobResult,
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
        return {
            handleJobResults: (payload) => {
                dispatch(setJobResults(payload));
            },
        };
    }
;


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(JobResultOptions));


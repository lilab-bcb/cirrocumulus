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
import {sortAndFilterJobResult, updateJob, updateTopNJobResult} from './JobResultsPanel';
import {
    createColorScale,
    INTERPOLATOR_SCALING_MIN_MAX_CATEGORY,
    INTERPOLATOR_SCALING_MIN_MAX_FEATURE,
    NATSORT
} from './util';

const styles = theme => ({
    formControl: {
        display: 'block',
        margin: theme.spacing(0, 1)
    },
    formControlInline: {
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
            min: '',
            max: ''
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
        jobResult.sizeScale = undefined;
        updateJob(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onSortChange = (event) => {
        const jobResult = this.getJobResult();
        jobResult.sortByGroup = event.target.value;
        updateTopNJobResult(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onMinUIChange = (value) => {
        this.setState({min: value});
    };
    onMaxUIChange = (value) => {
        this.setState({max: value});
    };

    onMinChange = (value) => {
        const jobResult = this.getJobResult();
        jobResult.options.min = value;
        jobResult.colorScale = undefined;
        updateJob(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onMaxChange = (value) => {
        const jobResult = this.getJobResult();
        jobResult.options.max = value;
        jobResult.colorScale = undefined;
        updateJob(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onOptions = (options) => {
        const jobResult = this.getJobResult();
        jobResult.options = Object.assign(jobResult.options, options);
        jobResult.colorScale = undefined;
        jobResult.sizeScale = undefined;
        this.handleJobResultDebounced(this.props.jobResultId);
    };

    onInterpolator = (value) => {
        const jobResult = this.getJobResult();
        const scale = jobResult.interpolator.scale;
        jobResult.interpolator = value;
        jobResult.interpolator.scale = scale;
        jobResult.colorScale = createColorScale(jobResult.interpolator).domain(jobResult.colorScale.domain()).clamp(true);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onColorScalingChange = (value) => {
        const jobResult = this.getJobResult();
        jobResult.colorScale = undefined;
        jobResult.interpolator.scale = value;
        updateJob(jobResult);
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onColorChanged = (event) => {
        const jobResult = this.getJobResult();
        jobResult.colorScale = undefined;
        jobResult.color = event.target.value;
        updateJob(this.getJobResult());
        this.props.handleJobResults(this.props.jobResults.slice());
    };

    onSizeChanged = (event) => {
        const jobResult = this.getJobResult();
        jobResult.sizeScale = undefined;
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
            by,
            colorScale,
            ntopUI,
            interpolator,
            sizeScale
        } = jobResult;

        const fields = jobResult.fields.slice();
        fields.sort(NATSORT);

        return <>
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
                    max={Math.min(100, jobResult.data.length)}
                    step={5}
                    style={{marginLeft: 10, width: '86%'}}
                    valueLabelDisplay="auto"
                    value={ntopUI}
                    onChange={this.onNtopChange}/>
                {jobResult.groups.length > 1 && <FormControl style={{display: 'block'}} className={classes.formControl}>
                    <InputLabel shrink={true}>Sort Features</InputLabel>
                    <Select
                        input={<Input size={"small"}/>}
                        onChange={this.onSortChange}
                        value={jobResult.sortByGroup}
                    >
                        {jobResult.columns.map(index => (
                            <MenuItem key={jobResult.groups[index]}
                                      value={jobResult.groups[index]}>{jobResult.groups[index]}</MenuItem>
                        ))}
                    </Select>
                </FormControl>}
            </div>
            <div style={{marginTop: 8}}></div>
            <Tooltip title="Filters are applied separately per cluster"><Typography
                component={"h2"}>Filters</Typography></Tooltip>

            {rowFilters.map(filter => {
                // [field, op, val, uiValue]
                const id = 'job_result' + filter[0];
                return <div key={filter[0]} style={{paddingTop: 8}}>
                    <FormControl className={classes.formControlInline}>
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

            <EditableColorScheme domain={colorScale.domain()}
                                 textColor={textColor}
                                 interpolator={interpolator}
                                 min={this.state.min}
                                 max={this.state.max}
                                 onMinChange={this.onMinChange}
                                 onMaxChange={this.onMaxChange}
                                 onMinUIChange={this.onMinUIChange}
                                 onMaxUIChange={this.onMaxUIChange}
                                 onInterpolator={this.onInterpolator}/>

            {jobResult.groups.length > 1 && <FormControl className={this.props.classes.formControl}>
                <InputLabel shrink={true}>Standardize</InputLabel>
                <Select
                    input={<Input size={"small"}/>}
                    onChange={event => this.onColorScalingChange(event.target.value)}
                    value={interpolator.scale}
                >
                    <MenuItem value={"none"} divider>(None)</MenuItem>
                    <MenuItem title={"Standardize features between 0 and 1"}
                              value={INTERPOLATOR_SCALING_MIN_MAX_FEATURE}>Feature</MenuItem>
                    <MenuItem title={"Standardize groups between 0 and 1"}
                              value={INTERPOLATOR_SCALING_MIN_MAX_CATEGORY}>Category</MenuItem>

                </Select>
            </FormControl>}

            <FormControl style={{marginTop: 16}} className={classes.formControl}>
                <InputLabel shrink={true}>Size</InputLabel>
                <Select
                    input={<Input size={"small"}/>}
                    onChange={this.onSizeChanged}
                    value={size}
                >
                    <MenuItem divider value={'none'}>{'(None)'}</MenuItem>
                    {fields.map(item => (
                        <MenuItem key={item} value={item}>{item}</MenuItem>
                    ))}
                </Select>
            </FormControl>
            <div style={{display: size === 'none' ? 'none' : '', marginTop: 6}}>
                <EditableSizeLegend sizeScale={sizeScale} textColor={textColor}
                                    onOptions={this.onOptions} showReversed={true}
                                    reversed={jobResult.sizeScaleReversed}
                                    onReversedChange={this.onSizeReversedChange}/>
            </div>
        </>;
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


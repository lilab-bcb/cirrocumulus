import {InputLabel, MenuItem, Select, Switch, Tooltip} from '@mui/material';
import FormControl from '@mui/material/FormControl';
import FormControlLabel from '@mui/material/FormControlLabel';
import Slider from '@mui/material/Slider';
import withStyles from '@mui/styles/withStyles';
import TextField from '@mui/material/TextField';
import Typography from '@mui/material/Typography';
import {debounce, findIndex} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setJobResults} from './actions';
import {EditableColorScheme} from './EditableColorScheme';
import {EditableSizeLegend} from './EditableSizeLegend';
import {sortAndFilterJobResult, updateJob, updateTopNJobResult} from './DotPlotJobResultsPanel';
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

class DotPlotJobResultOptions extends React.PureComponent {

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


    updateTopNJobResultDebounced = () => {
        const {jobResult} = this.props;
        updateTopNJobResult(jobResult);
        this.handleJobResults();
    };

    sortAndFilterDebounced = () => {
        const {jobResult} = this.props;
        sortAndFilterJobResult(jobResult);
        this.handleJobResults();
    };

    handleJobResultDebounced = () => {
        const {jobResult} = this.props;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onSizeReversedChange = (value) => {
        const {jobResult} = this.props;
        jobResult.sizeScaleReversed = value;
        jobResult.sizeScale = undefined;
        updateJob(jobResult);
        this.handleJobResults();
    };

    handleJobResults = () => {
        const {jobResult, jobResults} = this.props;
        const newJobResult = Object.assign({}, jobResult);
        const index = findIndex(jobResults, item => item.id === jobResult.id);
        if (index === -1) {
            throw new Error(jobResult.id + ' not found');
        }
        jobResults[index] = newJobResult; // TODO not redux style
        this.props.handleJobResults(jobResults.slice());
    };

    onSortChange = (event) => {
        const {jobResult} = this.props;
        jobResult.sortByGroup = event.target.value;
        updateTopNJobResult(jobResult);
        this.handleJobResults();
    };

    onMinUIChange = (value) => {
        this.setState({min: value});
    };
    onMaxUIChange = (value) => {
        this.setState({max: value});
    };

    onMinChange = (value) => {
        const {jobResult} = this.props;
        jobResult.options.min = value;
        jobResult.colorScale = undefined;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onMaxChange = (value) => {
        const {jobResult} = this.props;
        jobResult.options.max = value;
        jobResult.colorScale = undefined;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onOptions = (options) => {
        const {jobResult} = this.props;
        jobResult.options = Object.assign(jobResult.options, options);
        jobResult.colorScale = undefined;
        jobResult.sizeScale = undefined;
        this.handleJobResultDebounced(this.props.jobResultId);
    };

    onInterpolator = (value) => {
        const {jobResult} = this.props;
        const scale = jobResult.interpolator.scale;
        jobResult.interpolator = value;
        jobResult.interpolator.scale = scale;
        jobResult.colorScale = createColorScale(jobResult.interpolator).domain(jobResult.colorScale.domain()).clamp(true);
        this.handleJobResults();
    };

    onColorScalingChange = (value) => {
        const {jobResult} = this.props;
        jobResult.colorScale = undefined;
        jobResult.interpolator.scale = value;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onColorChanged = (event) => {
        const {jobResult} = this.props;
        jobResult.colorScale = undefined;
        jobResult.color = event.target.value;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onSizeChanged = (event) => {
        const {jobResult} = this.props;
        jobResult.sizeScale = undefined;
        jobResult.size = event.target.value;
        updateJob(jobResult);
        this.handleJobResults();
    };

    onNtopChange = (event, value) => {
        const {jobResult} = this.props;
        jobResult.ntopUI = value;
        jobResult.ntop = value;
        this.setState((prevState, props) => ({
            forceUpdate: !prevState.forceUpdate
        }));
        this.updateTopNJobResultDebounced();
    };

    onByChanged = (event) => {
        const {jobResult} = this.props;
        jobResult.by = event.target.value;
        sortAndFilterJobResult(jobResult);
        this.handleJobResults();
    };

    onByAscendingChange = (event) => {
        const {jobResult} = this.props;
        jobResult.byAscending = event.target.checked;
        sortAndFilterJobResult(jobResult);
        this.handleJobResults();
    };

    onOperationChanged = (filter, value) => {
        const {jobResult} = this.props;
        filter[1] = value;
        sortAndFilterJobResult(jobResult);
        this.handleJobResults();
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
            jobResult,
            classes,
            textColor
        } = this.props;


        const {
            rowFilters,
            color,
            size,
            by,
            byAscending,
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
                    <InputLabel>By</InputLabel>
                    <Select
                        label={"By"}
                        size={"small"}
                        onChange={this.onByChanged}
                        value={by}
                    >
                        {fields.map(item => (
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl>
                <div><FormControlLabel
                    control={
                        <Switch
                            value={"byAscending"}
                            checked={byAscending}
                            onChange={this.onByAscendingChange}
                        />
                    }
                    label="Ascending"
                /></div>
                <InputLabel style={{marginLeft: 8, marginTop: 8}}>Number of Features</InputLabel>
                <Slider
                    min={5}
                    max={Math.min(100, jobResult.data.length)}
                    step={5}
                    style={{marginLeft: 10, width: '86%'}}
                    valueLabelDisplay="auto"
                    value={ntopUI}
                    onChange={this.onNtopChange}/>
                {jobResult.groups.length > 1 && <FormControl style={{display: 'block'}} className={classes.formControl}>
                    <InputLabel>Sort Features</InputLabel>
                    <Select
                        label={"Sort Features"}
                        size={"small"}
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
                            label={filter[0]}
                            size={"small"}
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
                    <TextField size={"small"}
                               onChange={event => this.onValueChange(filter, event.target.value)} value={filter[3]}
                               style={{maxWidth: 60, verticalAlign: 'bottom'}}/>
                </div>;
            })}


            <Typography style={{marginTop: 16}}
                        component={"h2"}>View Options</Typography>

            <FormControl className={classes.formControl}>
                <InputLabel>Color</InputLabel>
                <Select
                    label={"Color"}
                    size={"small"}
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
                <InputLabel>Standardize</InputLabel>
                <Select
                    label={"Standardize"}
                    size={"small"}
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
                <InputLabel>Size</InputLabel>
                <Select
                    label={"Size"}
                    size={"small"}
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
            jobResults: state.jobResults
        };
    }
;
const mapDispatchToProps = (dispatch) => {
        return {
            handleJobResults: (payload) => {
                dispatch(setJobResults(payload));
            }
        };
    }
;


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(DotPlotJobResultOptions));


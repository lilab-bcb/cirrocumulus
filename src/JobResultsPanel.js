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
    TableRow,
    Tooltip
} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import Divider from '@material-ui/core/Divider';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import withStyles from '@material-ui/core/styles/withStyles';
import Typography from '@material-ui/core/Typography';
import CloudDownloadIcon from '@material-ui/icons/CloudDownload';
import DeleteIcon from '@material-ui/icons/Delete';
import {scaleLinear} from 'd3-scale';
import {find} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {deleteJobResult, setJobResult, setSearchTokensDirectly, setTab} from './actions';
import {createFilterFunction} from './dataset_filter';
import {intFormat} from './formatters';
import {
    createColorScale,
    FEATURE_TYPE,
    getInterpolator,
    INTERPOLATOR_SCALING_MIN_MAX_CATEGORY,
    INTERPOLATOR_SCALING_MIN_MAX_FEATURE,
    INTERPOLATOR_SCALING_NONE,
    NATSORT,
    scaleConstantRange
} from './util';

import DotPlotTable from './DotPlotTable';
import {COMPARE_ACTIONS} from './job_config';

const styles = theme => ({
    toolbar: {
        '& hr': {
            margin: theme.spacing(0, 0.5)
        }
    }
});


const DEFAULT_DE_INTERPOLATOR = 'RdBu';


export function updateJob(jobResult) {
    if (jobResult.options === undefined) {
        jobResult.options = {};
    }
    if (jobResult.interpolator === undefined) {
        jobResult.interpolator = {
            name: DEFAULT_DE_INTERPOLATOR,
            value: getInterpolator(DEFAULT_DE_INTERPOLATOR),
            reversed: true,
            scale: INTERPOLATOR_SCALING_NONE
        };
    }
    const groups = jobResult.groups;
    const data = jobResult.data;

    if (jobResult.columns === undefined) {
        const indices = new Array(groups.length);
        for (let i = 0, n = jobResult.groups.length; i < n; i++) {
            indices[i] = i;
        }
        jobResult.columns = indices;
        jobResult.columns.sort((a, b) => {
            return NATSORT(jobResult.groups[a], jobResult.groups[b]);
        });
    }

    if (jobResult.rowFilters === undefined) {
        const filters = [];
        jobResult.fields.forEach(field => {
            const fieldLowercase = field.toLowerCase();
            filters.push([field, fieldLowercase.indexOf('pval') !== -1 || fieldLowercase.indexOf('p_val') !== -1 || fieldLowercase.indexOf('qval') !== -1 || fieldLowercase.indexOf('q_val') !== -1 || field.indexOf('fdr') !== -1 ? '<' : '>', NaN, '']);
        });

        jobResult.rowFilters = filters;
    }
    if (jobResult.ntop === undefined) {
        const ntop = Math.min(10, jobResult.data.length);
        jobResult.ntop = ntop;
        jobResult.ntopUI = ntop;
    }
    if (jobResult.sortedRows === undefined) {
        const indices = new Array(data.length);
        for (let i = 0, n = data.length; i < n; i++) {
            indices[i] = i;
        }
        jobResult.sortedRows = [];
        for (let i = 0, n = groups.length; i < n; i++) {
            jobResult.sortedRows.push(indices);
        }
    }
    if (jobResult.byAscending === undefined) {
        jobResult.byAscending = false;
    }
    if (jobResult.by === undefined) {
        jobResult.by = jobResult.fields[0];
        for (let i = 0; i < jobResult.fields.length; i++) {
            if (jobResult.fields[i].toLowerCase().indexOf('score') !== -1) {
                jobResult.by = jobResult.fields[i];
                break;
            }
        }
    }
    if (jobResult.sortByGroup === undefined) {
        jobResult.sortByGroup = jobResult.groups[0];
    }
    if (jobResult.color === undefined) {
        jobResult.color = jobResult.fields[0];
    }
    if (jobResult.size === undefined) {
        jobResult.size = jobResult.fields[0];
    }
    if (jobResult.rowSortOrder === undefined) {
        jobResult.rowSortOrder = [];
    }
    if (jobResult.rows === undefined) {
        sortAndFilterJobResult(jobResult);
    }


    function getRange(field) {
        let min = Number.MAX_VALUE;
        let max = -Number.MAX_VALUE;
        for (let i = 0, n = data.length; i < n; i++) {
            for (let j = 0; j < groups.length; j++) {
                const group = groups[j];
                const fullField = group + ':' + field;
                const value = data[i][fullField];
                if (value != null && !isNaN(value)) {
                    min = Math.min(min, value);
                    max = Math.max(max, value);
                }
            }
        }
        return [min, max];
    }

    // color='logfoldchanges', size='pvals_adj',
    if (jobResult.colorScale === undefined) {
        let domain;
        if (jobResult.interpolator.scale !== INTERPOLATOR_SCALING_NONE) {
            domain = [0, 1];
        } else {
            domain = getRange(jobResult.color);
        }
        // if (isNaN(jobResult.options.min) && isNaN(jobResult.options.max) && domain[0] < 0 && domain[1] > 0) {
        //     const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
        //     domain[0] = -max;
        //     domain[1] = max;
        // }
        if (!isNaN(jobResult.options.min)) {
            domain[0] = jobResult.options.min;
        }
        if (!isNaN(jobResult.options.max)) {
            domain[1] = jobResult.options.max;
        }
        jobResult.colorScale = createColorScale(jobResult.interpolator).domain(domain);
    }
    if (jobResult.sizeScaleReversed === undefined) {
        jobResult.sizeScaleReversed = jobResult.size != null && jobResult.size.indexOf('pval') !== -1;
    }
    if (jobResult.sizeScale === undefined) {
        if (jobResult.size !== 'none') {
            let domain = getRange(jobResult.size);
            if (!isNaN(jobResult.options.minSize)) {
                domain[0] = jobResult.options.minSize;
            }
            if (!isNaN(jobResult.options.maxSize)) {
                domain[1] = jobResult.options.maxSize;
            }
            jobResult.sizeScale = scaleLinear().domain(domain).range(jobResult.sizeScaleReversed ? [18, 2] : [2, 18]).clamp(true);
        } else {
            jobResult.sizeScale = scaleConstantRange(18);
        }
    }
}

export function sortAndFilterJobResult(jobResult) {
    const by = jobResult.by;
    const byAscending = jobResult.byAscending;
    const groups = jobResult.groups;
    const ngroups = groups.length;
    jobResult.sortedRows = [];
    // sort each group
    for (let groupIndex = 0; groupIndex < ngroups; groupIndex++) {
        const indices = [];
        for (let i = 0, n = jobResult.data.length; i < n; i++) {
            indices.push(i);
        }

        const group = groups[groupIndex];
        const byField = group + ':' + by;
        indices.sort((a, b) => {
            const val1 = jobResult.data[a][byField];
            const val2 = jobResult.data[b][byField];
            const aNaN = (val1 == null || isNaN(val1));
            const bNaN = (val2 == null || isNaN(val2));
            if (aNaN && bNaN) {
                return 0;
            }
            if (byAscending) {
                if (aNaN) {
                    return -1;
                }
                if (bNaN) {
                    return 1;
                }
                return val1 - val2;
            } else {
                if (aNaN) {
                    return 1;
                }
                if (bNaN) {
                    return -1;
                }
                return val2 - val1;
            }
        });
        jobResult.sortedRows.push(indices);
    }

    const rowFilters = jobResult.rowFilters;
    const rowFilterFunctions = [];
    rowFilters.forEach(filter => {
        if (!isNaN(filter[2])) {
            rowFilterFunctions.push({field: filter[0], f: createFilterFunction(filter)});
        }
    });

    // apply filters separately per cluster
    const nFilters = rowFilterFunctions.length;
    // filter
    jobResult.sortedFilteredRows = [];
    for (let groupIndex = 0; groupIndex < ngroups; groupIndex++) {
        const group = groups[groupIndex];
        const filteredIndices = [];
        const sortedIndices = jobResult.sortedRows[groupIndex];
        for (let i = 0, n = sortedIndices.length; i < n; i++) {
            let passesFilter = true;
            const rowIndex = sortedIndices[i];
            for (let filterIndex = 0; filterIndex < nFilters; filterIndex++) {
                const filter = rowFilterFunctions[filterIndex];
                const field = group + ':' + filter.field;
                const value = jobResult.data[rowIndex][field];
                if (!filter.f(value) || value == null || isNaN(value)) {
                    passesFilter = false;
                    break;
                }
            }
            if (passesFilter) {
                filteredIndices.push(rowIndex);
            }
        }
        jobResult.sortedFilteredRows.push(filteredIndices);
    }
    updateTopNJobResult(jobResult);
}


export function updateTopNJobResult(jobResult) {
    let indices = new Set();
    const index = jobResult.groups.indexOf(jobResult.sortByGroup);
    const groupOrder = [index];
    for (let i = 0; i < jobResult.groups.length; i++) {
        if (i !== index) {
            groupOrder.push(i);
        }
    }
    for (let i = 0; i < jobResult.groups.length; i++) {
        const groupIndices = jobResult.sortedFilteredRows[groupOrder[i]];
        const ntop = Math.min(jobResult.ntop, groupIndices.length);
        for (let featureIndex = 0; featureIndex < ntop; featureIndex++) {
            indices.add(groupIndices[featureIndex]);
        }
    }
    jobResult.rows = Array.from(indices);
}

class JobResultsPanel extends React.PureComponent {


    constructor(props) {
        super(props);
        this.state = {showDialog: false};
    }

    onMouseMove = (event) => {
        this.props.setTooltip(event.target.dataset.title);
    };

    onMouseOut = (event) => {
        this.props.setTooltip('');
    };

    onBrowseJobs = (event) => {
        this.setState({showDialog: true});
    };

    onCloseJobs = () => {
        this.setState({showDialog: false});
    };

    onSelectJob = (id) => {
        this.props.setJobResult(id);
        this.setState({showDialog: false});
    };

    onDeleteJob = (event, job) => {
        event.stopPropagation();
        this.setState({browseJob: job});
    };

    onDeleteJobOK = (id) => {
        this.props.onDeleteJob(id);
        this.setState({browseJob: null});
    };

    onDeleteJobCancel = () => {
        this.setState({browseJob: null});
    };
    getJobResult = (id) => {
        const {jobResults, jobResultId} = this.props;
        if (id == null) {
            id = jobResultId;
        }
        return find(jobResults, item => item.id === id);
    };

    exportJobResult = (event) => {
        const jobResult = this.getJobResult();
        const output = [];
        output.push('id');
        jobResult.columns.forEach(columnIndex => {
            const group = jobResult.groups[columnIndex];
            jobResult.fields.forEach(field => {
                output.push('\t');
                output.push(group + ':' + field);
            });
        });
        output.push('\n');
        for (let i = 0; i < jobResult.data.length; i++) {
            output.push(jobResult.data[i].index);
            jobResult.columns.forEach(columnIndex => {
                const group = jobResult.groups[columnIndex];
                jobResult.fields.forEach(field => {
                    const fieldName = group + ':' + field;
                    const value = jobResult.data[i][fieldName];
                    output.push('\t');
                    output.push(value);
                });
            });
            output.push('\n');
        }
        const blob = new Blob([output.join('')], {
            type: 'text/plain;charset=utf-8'
        });
        window.saveAs(blob, jobResult.name + '.tsv');
    };


    toggleAll = (event, isSelected) => {
        event.stopPropagation();
        isSelected = !isSelected;

        const jobResult = this.getJobResult();
        const features = new Set();
        for (let i = 0; i < jobResult.rows.length; i++) {
            const feature = jobResult.data[jobResult.rows[i]]['index'];
            features.add(feature);
        }

        let searchTokens = this.props.searchTokens;
        if (!isSelected) {
            searchTokens = searchTokens.filter(token => !features.has(token.value));
        } else {
            const dataset = this.props.dataset;
            features.forEach(feature => {
                let found = false;
                for (let i = 0; i < searchTokens.length; i++) {
                    if (searchTokens[i].value === feature) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    searchTokens.push({
                        value: feature,
                        type: dataset.obs.indexOf(feature) !== -1 ? FEATURE_TYPE.OBS : FEATURE_TYPE.X
                    });
                }
            });
        }
        this.props.onSearchTokens(searchTokens.slice());
    };

    onRowClick = (event, row) => {
        event.stopPropagation();
        let searchTokens = this.props.searchTokens;
        const jobResult = this.getJobResult();
        const feature = jobResult.data[row]['index'];
        let index = -1;
        for (let i = 0; i < searchTokens.length; i++) {
            if (searchTokens[i].value === feature) {
                index = i;
                break;
            }
        }
        if (index === -1) {
            const dataset = this.props.dataset;
            searchTokens.push({
                value: feature,
                type: dataset.obs.indexOf(feature) !== -1 ? FEATURE_TYPE.OBS : FEATURE_TYPE.X
            });
        } else {
            searchTokens.splice(index, 1);
        }
        this.props.onSearchTokens(searchTokens.slice());
    };

    sortColumns = () => {
        const jobResult = this.getJobResult();
        const sortOrder = jobResult.columnSortOrder;
        const groups = jobResult.groups;
        const data = jobResult.data;
        const indices = new Array(groups.length);
        for (let i = 0, n = groups.length; i < n; i++) {
            indices[i] = i;
        }
        const nsortFields = sortOrder.length;
        if (nsortFields.length > 0) { // sort by gene
            const sortFieldRows = [];
            for (let i = 0; i < nsortFields; i++) {
                const field = sortOrder[i].field;
                let index;
                for (let j = 0; j < data.length; j++) {
                    if (data[j]['index'] === field) {
                        index = j;
                        break;
                    }
                }
                sortFieldRows.push(index);
            }

            indices.sort((a, b) => {
                for (let i = 0; i < nsortFields; i++) {
                    const row = sortFieldRows[i];
                    const val = data[row][jobResult.color] - data[row][jobResult.color];
                    if (val !== 0) {
                        return val;
                    }
                    const val2 = data[row][jobResult.size] - data[row][jobResult.size];
                    if (val2 !== 0) {
                        return val2;
                    }
                }
                return a - b;
            });
        }
        jobResult.columns = indices;
    };

    render() {
        const {email, jobResultId, jobResults, classes, searchTokens} = this.props;
        const jobResult = this.getJobResult();
        const selectedFeatures = new Set();
        searchTokens.forEach(token => {
            if (token.type === FEATURE_TYPE.X) {
                selectedFeatures.add(token.value);
            }
        });

        let data;
        let color;
        let by;
        let size;
        let groups;
        let rotateHeaders = false;
        let rows;
        let valueScale = null;
        let isSizeScaled = true;
        let domains;
        let tooltipFields = null;
        let headerHeight = 0;
        let headerWidth = 0;
        // let isSingleComparison = false;
        if (jobResult != null) {
            const columns = jobResult.columns;
            // isSingleComparison = columns.length === 1;
            isSizeScaled = jobResult.size !== 'none';
            data = jobResult.data;
            color = jobResult.color;
            by = jobResult.by;
            size = jobResult.size;
            groups = jobResult.groups;

            // if (isSingleComparison) {
            //     // min/max for each field
            //     jobResult.fields.forEach(field => {
            //         let min = Number.MAX_VALUE;
            //         let max = -Number.MAX_VALUE;
            //         const fieldName = 'comparison:' + field;
            //         rows.forEach(row => {
            //             const value = data[row][fieldName];
            //             if (value != null && !isNaN(value)) {
            //                 min = value < min ? value : min;
            //                 max = value > max ? value : max;
            //             }
            //         });
            //         domains.push([min, max]);
            //     });
            // }
            tooltipFields = [];
            tooltipFields.push(by);
            rows = jobResult.rows;
            if (color !== by) {
                tooltipFields.push(color);
            }
            if (size !== by && isSizeScaled && size !== color) {
                tooltipFields.push(size);
            }

            jobResult.fields.forEach(field => {
                if (tooltipFields.indexOf(field) === -1) {
                    tooltipFields.push(field);
                }
            });


            for (let i = 0; i < groups.length; i++) {
                if (groups[i].length > 2) {
                    rotateHeaders = true;
                    break;
                }
            }

            if (rotateHeaders) {
                const d = document.createElement('span');
                d.style.position = 'absolute';
                d.style.left = '-1000000px';
                d.className = '.MuiTableCell-head .MuiTableCell-root';
                document.body.append(d);
                for (let i = 0; i < groups.length; i++) {
                    d.innerText = groups[i];
                    headerWidth = Math.max(headerWidth, 2 + d.getBoundingClientRect().width);
                }
                d.remove();
                headerWidth += 2; // prevent overflow
                headerWidth = Math.min(headerWidth, 300);
                headerHeight = Math.cos(45) * headerWidth;
            }
            if (jobResult.interpolator.scale !== INTERPOLATOR_SCALING_NONE) {
                valueScale = scaleLinear().range([0, 1]);
                domains = [];
                // can scale colors globally, by row, or by column
                if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                    rows.forEach(row => {
                        let min = Number.MAX_VALUE;
                        let max = -Number.MAX_VALUE;
                        columns.forEach(column => {
                            const group = groups[column];
                            const colorField = group + ':' + color;
                            const colorValue = data[row][colorField];
                            if (colorValue != null && !isNaN(colorValue)) {
                                min = colorValue < min ? colorValue : min;
                                max = colorValue > max ? colorValue : max;
                            }
                        });
                        domains.push([min, max]);
                    });
                } else if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
                    columns.forEach(column => {
                        let min = Number.MAX_VALUE;
                        let max = -Number.MAX_VALUE;
                        rows.forEach(row => {
                            const group = groups[column];
                            const colorField = group + ':' + color;
                            const colorValue = data[row][colorField];
                            if (colorValue != null && !isNaN(colorValue)) {
                                min = colorValue < min ? colorValue : min;
                                max = colorValue > max ? colorValue : max;
                            }
                        });
                        domains.push([min, max]);
                    });
                }
            }

        }

        let showJobStatus = false;
        for (let i = 0; i < jobResults.length; i++) {
            const isPrecomputed = ('' + jobResults[i].id).startsWith('cirro-');
            if (!isPrecomputed) {
                showJobStatus = true;
                break;
            }
        }
        const showBrowseJobs = (jobResultId == null && jobResults.length > 0) || jobResults.length > 1;
        const columns = jobResult ? jobResult.columns.map(column => groups[column]) : null;
        const isRowSelected = (row) => selectedFeatures.has(data[row].index);
        const getRowId = (row) => data[row].index;
        const getTooltip = (row, column) => {
            let title = [];
            tooltipFields.forEach(field => {
                const value = data[row][column + ':' + field];
                title.push(field + ': ' + value);
            });
            title = title.join(', ');
            return title;
        };


        const getColor = (row, column) => data[row][column + ':' + color];
        const getSize = (row, column) => data[row][column + ':' + size];

        const rowStart = (row, rowIndex) => {
            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                valueScale.domain(domains[rowIndex]);
            }
        };
        const columnStart = (column, columnIndex) => {
            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
                valueScale.domain(domains[columnIndex]);
            }
        };
        const jobTypeToName = {};
        COMPARE_ACTIONS.forEach(action => jobTypeToName[action.jobType] = action.title);
        return <>
            <Box color="text.primary">
                <Grid container alignItems="center" className={classes.toolbar}>
                    {showBrowseJobs &&
                    <Button size={"small"} onClick={this.onBrowseJobs} variant="outlined"
                            color="primary">Browse All Results</Button>}
                    {showBrowseJobs && jobResult && <Divider orientation="vertical" flexItem/>}
                    {jobResult && <Tooltip title={"Export"}><IconButton edge={false} size={'small'} aria-label="Export"
                                                                        onClick={this.exportJobResult}><CloudDownloadIcon/></IconButton></Tooltip>}
                </Grid>
                {jobResult && <>
                    <Typography
                        style={{marginBottom: headerHeight + 8}}
                        component={"h3"}
                        color={"textPrimary"}><b>{jobResult.name}</b>&nbsp;
                        <small>{intFormat(rows.length) + ' / ' + intFormat(jobResult.data.length) + ' features'}</small></Typography>
                    <div style={{paddingTop: 6}}>
                        {rows.length > 0 &&
                        <DotPlotTable isRowSelected={isRowSelected}
                                      rows={rows}
                                      columns={columns}
                                      onMouseMove={this.onMouseMove}
                                      onMouseOut={this.onMouseOut}
                                      getTooltip={getTooltip}
                                      sizeScale={jobResult.sizeScale}
                                      valueScale={valueScale}
                                      rowStart={rowStart}
                                      getColor={getColor}
                                      getSize={getSize}
                                      getRowId={getRowId}
                                      columnStart={columnStart}
                                      onRowClick={this.onRowClick}
                                      colorScale={jobResult.colorScale}
                                      toggleAll={this.toggleAll}
                                      rotateHeaders={rotateHeaders}
                                      headerHeight={headerHeight}
                                      headerWidth={headerWidth}/>}
                    </div>
                </>
                }
            </Box>
            <Dialog
                open={this.state.browseJob != null}
                onClose={this.onDeleteJobCancel}
                aria-labelledby="confirm-dialog-title"
                aria-describedby="conform-dialog-description"
            >
                <DialogTitle id="confirm-dialog-title">Confirm Delete</DialogTitle>
                <DialogContent>
                    <DialogContentText id="confirm-dialog-description">
                        Are you sure you want to delete {this.state.browseJob ? this.state.browseJob.name : ''}?
                    </DialogContentText>
                </DialogContent>
                <DialogActions>
                    <Button onClick={this.onDeleteJobCancel}>
                        Cancel
                    </Button>
                    <Button variant="contained" onClick={e => this.onDeleteJobOK(this.state.browseJob.id)}
                            color="primary"
                            autoFocus>
                        OK
                    </Button>
                </DialogActions>
            </Dialog>
            <Dialog onClose={this.onCloseJobs} aria-labelledby="job-results-title"
                    open={this.state.showDialog && jobResults.length > 0}>
                <DialogTitle id="job-results-title" onClose={this.onCloseJobs}>
                    Results
                </DialogTitle>
                <DialogContent>
                    <Table size="small" stickyHeader={true}>
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
                                const isComplete = isPrecomputed || jobResult.status === 'complete';

                                const status = isPrecomputed ? 'complete' : jobResult.status;
                                const isJobOwner = email == jobResult.email || (email === null && jobResult.email === '');
                                const jobType = jobTypeToName[jobResult.type];
                                // const date = isPrecomputed ? '' : jobResult.submitted;
                                return <TableRow key={jobResult.id}
                                                 className={classes.deleteTr}
                                                 hover
                                                 selected={jobResult.id === jobResultId}
                                                 onClick={isComplete ? (event) => this.onSelectJob(jobResult.id) : null}
                                                 role="checkbox"
                                                 tabIndex={-1}>
                                    <TableCell>{text}</TableCell>
                                    <TableCell>{jobType}</TableCell>
                                    {showJobStatus && <TableCell>{status}
                                        {isJobOwner && !isPrecomputed &&
                                        <IconButton edge="end" aria-label="delete"
                                                    onClick={(event) => this.onDeleteJob(event, jobResult)}>
                                            <DeleteIcon/>
                                        </IconButton>}</TableCell>}
                                </TableRow>;
                            })}
                        </TableBody>
                    </Table>
                </DialogContent>
            </Dialog>
        </>;
    };
}


const mapStateToProps = state => {
        return {
            dataset: state.dataset,
            email: state.email,
            jobResultId: state.jobResult,
            jobResults: state.jobResults,
            searchTokens: state.searchTokens,
            tab: state.tab
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
        return {
            onSearchTokens: (payload) => {
                dispatch(setSearchTokensDirectly(payload));
            },
            onDeleteJob: (payload) => {
                dispatch(deleteJobResult(payload));
            },
            setTab: (payload) => {
                dispatch(setTab(payload));
            },
            setJobResult: (payload) => {
                dispatch(setJobResult(payload));
            }
        };
    }
;


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(JobResultsPanel));

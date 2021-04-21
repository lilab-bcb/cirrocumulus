import {
    Checkbox,
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
import {intFormat, numberFormat2f} from './formatters';
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

const styles = theme => ({
    dot: {
        borderRadius: '50%',
        display: 'inlineBlock',
        border: '1px solid lightgray'
    },
    toolbar: {
        '& hr': {
            margin: theme.spacing(0, 0.5),
        }
    },
    table: {
        borderCollapse: 'collapse',
        position: 'relative',
        width: 'unset'
    },
    tr: {
        cursor: 'pointer',
    },
    deleteTr: {
        cursor: 'pointer',
        '& span': {
            display: 'none',
        },
        '&:hover span': {
            display: 'block',
            position: 'absolute',
            right: 0,
            top: 0
        }
    },
    rotateHeader: {
        top: 0,
        padding: 1,
        background: 'transparent',
        whiteSpace: 'nowrap'
    },
    checkbox: {
        padding: 0
    },
    rotateHeaderDiv: {
        /* place div at bottom left of the th parent */
        position: 'absolute',
        bottom: 0,
        left: 0,
        /* Make sure short labels still meet the corner of the parent otherwise you'll get a gap */
        textAlign: 'left',
        /* Move the top left corner of the span's bottom-border to line up with the top left corner of the td's border-right border so that the border corners are matched
         * Rotate 315 (-45) degrees about matched border corners */
        transform: 'translate(calc(50%),0) rotate(315deg)',
        transformOrigin: '0% calc(50%)',
        width: '100%',
    },
    rotateHeaderSpan: {
        position: 'absolute',
        bottom: 0,
        left: 0,
        pointerEvents: 'none',
        overflow: 'hidden',
        textOverflow: 'ellipsis'
    },
    td: {
        padding: 1,
        whiteSpace: 'nowrap',
        width: 27,
        minWidth: 27,
        maxWidth: 27
    },
    rowHeader: {
        padding: 1,
        whiteSpace: 'nowrap',
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
            filters.push([field, '>', NaN, '']);
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
            if (aNaN) {
                return 1;
            }
            if (bNaN) {
                return -1;
            }
            return val2 - val1;
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
    handleClick = (event, row) => {
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
        // const name = jobResult.name;
        // const params = jobResult.params;
        let data;
        let color;
        let by;
        let size;
        let groups;
        let rotateHeaders = false;
        let rows;
        let headerWidth = 0;
        let headerHeight = 0;
        let maxSize;
        let valueScale = null;
        let isSizeScaled = true;
        let domains;
        let tooltipFields = null;
        let selectAllChecked = true;
        if (jobResult != null) {
            // const name = jobResult.name;
            // const params = jobResult.params;
            isSizeScaled = jobResult.size !== 'none';
            data = jobResult.data;
            color = jobResult.color;
            by = jobResult.by;
            size = jobResult.size;
            groups = jobResult.groups;
            tooltipFields = [];
            tooltipFields.push(by);

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

            if (jobResult.interpolator.scale !== INTERPOLATOR_SCALING_NONE) {
                valueScale = scaleLinear().range([0, 1]);
                domains = [];
                if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                    jobResult.rows.forEach(row => {
                        let min = Number.MAX_VALUE;
                        let max = -Number.MAX_VALUE;
                        jobResult.columns.forEach(column => {
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
                    jobResult.columns.forEach(column => {
                        let min = Number.MAX_VALUE;
                        let max = -Number.MAX_VALUE;
                        jobResult.rows.forEach(row => {
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

            for (let i = 0; i < groups.length; i++) {
                if (groups[i].length > 2) {
                    rotateHeaders = true;
                    break;
                }
            }
            rows = jobResult.rows;

            for (let i = 0; i < rows.length; i++) {
                const id = data[rows[i]]['index'];
                const selected = selectedFeatures.has(id);
                if (!selected) {
                    selectAllChecked = false;
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
                    headerWidth = Math.max(headerWidth, d.getBoundingClientRect().width);
                }
                d.remove();
                headerWidth += 2; // prevent overflow
                headerWidth = Math.min(headerWidth, 300);
                headerHeight = Math.cos(45) * headerWidth;
            }
            maxSize = Math.max(jobResult.sizeScale.range()[0], jobResult.sizeScale.range()[1]);
        }

        function formatNumber(val) {
            return (val == null || isNaN(val)) ? "NaN" : numberFormat2f(val);
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
                        {jobResult.rows.length > 0 &&
                        <Table onMouseMove={this.onMouseMove} onMouseOut={this.onMouseOut} stickyHeader={true}
                               className={classes.table}>
                            <TableHead>
                                <TableRow>
                                    <TableCell className={classes.rowHeader} component="th"
                                               style={{backgroundColor: 'unset', textAlign: 'left'}}
                                               key={'__id'}><Checkbox
                                        onClick={(event) => this.toggleAll(event, selectAllChecked)}
                                        className={classes.checkbox}
                                        checked={selectAllChecked}/></TableCell>
                                    {jobResult.columns.map(column => {
                                        const group = groups[column];
                                        if (rotateHeaders) {
                                            return <TableCell
                                                title={group}
                                                className={classes.rotateHeader}
                                                key={group}>
                                                <div className={classes.rotateHeaderDiv}><span
                                                    style={{width: headerWidth}}
                                                    className={classes.rotateHeaderSpan}>{group}</span></div>
                                            </TableCell>;
                                        } else {
                                            return <TableCell style={{backgroundColor: 'unset'}}
                                                              className={classes.td}
                                                              key={group}>{group}</TableCell>;
                                        }
                                    })}
                                </TableRow>
                            </TableHead>
                            <TableBody>
                                {rows.map((row, rowIndex) => {
                                    const id = data[row]['index'];
                                    const selected = selectedFeatures.has(id);
                                    if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                                        valueScale.domain(domains[rowIndex]);
                                    }
                                    return <TableRow
                                        className={classes.tr}
                                        hover
                                        onClick={(event) => this.handleClick(event, row)}
                                        role="checkbox"
                                        tabIndex={-1}
                                        key={row}
                                    ><TableCell className={classes.rowHeader} component="th"
                                                key={'id'}><Checkbox className={classes.checkbox}
                                                                     checked={selected}/>{id}</TableCell>
                                        {jobResult.columns.map((column, columnIndex) => {
                                            const group = groups[column];
                                            const colorField = group + ':' + color;
                                            const sizeField = group + ':' + size;
                                            const colorValue = data[row][colorField];
                                            const sizeValue = data[row][sizeField];
                                            let title = [];
                                            tooltipFields.forEach(field => {
                                                const value = data[row][group + ':' + field];
                                                title.push(field + ': ' + value);
                                            });
                                            title = title.join(', ');
                                            if (colorValue == null || isNaN(colorValue)) {
                                                return <TableCell className={classes.td} data-title={title}
                                                                  key={group}/>;
                                            }
                                            const diameter = jobResult.sizeScale(sizeValue);
                                            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
                                                valueScale.domain(domains[columnIndex]);
                                            }
                                            const colorValueScaled = valueScale ? valueScale(colorValue) : colorValue;
                                            const backgroundColor = jobResult.colorScale(colorValueScaled);
                                            return <TableCell className={classes.td} data-title={title}
                                                              key={group}>
                                                <div className={classes.dot}
                                                     style={{
                                                         pointerEvents: 'none',
                                                         marginLeft: (maxSize - diameter) / 2,
                                                         width: diameter,
                                                         height: diameter,
                                                         backgroundColor: backgroundColor
                                                     }}></div>
                                            </TableCell>;
                                        })}
                                    </TableRow>;
                                })}
                            </TableBody>
                        </Table>}
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
                    open={this.state.showDialog}>
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
                                const jobType = jobResult.type === 'de' ? 'Differential Expression' : 'Correlation';
                                const status = isPrecomputed ? 'complete' : jobResult.status;
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
                                        {email == jobResult.email && !isPrecomputed &&
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
            },
        };
    }
;


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(JobResultsPanel));

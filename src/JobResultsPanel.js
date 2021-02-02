import {
    Checkbox,
    Dialog,
    DialogContent,
    DialogTitle,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow
} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import withStyles from '@material-ui/core/styles/withStyles';
import Typography from '@material-ui/core/Typography';
import {scaleLinear} from 'd3-scale';
import {find} from 'lodash';
import natsort from 'natsort';
import React from 'react';
import {connect} from 'react-redux';
import {setJobResult, setSearchTokensDirectly, setTab} from './actions';
import {createFilterFunction} from './dataset_filter';
import {intFormat, numberFormat2f} from './formatters';
import {createColorScale, getInterpolator, X_SEARCH_TOKEN} from './util';

const styles = theme => ({
    dot: {
        borderRadius: '50%',
        display: 'inlineBlock',
        border: '1px solid lightgray'
    },
    visuallyHidden: {
        border: 0,
        clip: 'rect(0 0 0 0)',
        height: 1,
        margin: -1,
        overflow: 'hidden',
        padding: 0,
        position: 'absolute',
        top: 20,
        width: 1,
    },
    table: {
        borderCollapse: 'collapse',
        position: 'relative',
        width: 'unset'
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
    if (jobResult.options == null) {
        jobResult.options = {};
    }
    if (jobResult.interpolator == null) {
        jobResult.interpolator = {
            name: DEFAULT_DE_INTERPOLATOR,
            value: getInterpolator(DEFAULT_DE_INTERPOLATOR),
            reversed: true
        };
    }
    const groups = jobResult.groups;
    const data = jobResult.data;

    if (jobResult.columns == null) {
        const sorter = natsort({insensitive: true});
        const indices = new Array(groups.length);
        for (let i = 0, n = jobResult.groups.length; i < n; i++) {
            indices[i] = i;
        }
        jobResult.columns = indices;
        jobResult.columns.sort((a, b) => {
            return sorter(jobResult.groups[a], jobResult.groups[b]);
        });
    }

    if (jobResult.rowFilters == null) {
        const filters = [];
        jobResult.fields.forEach(field => {
            filters.push([field, '>', NaN, '']);
        });

        jobResult.rowFilters = filters;
    }
    if (jobResult.ntop == null) {
        jobResult.ntop = 10;
        jobResult.ntopUI = 10;
    }
    if (jobResult.sortedRows == null) {
        const indices = new Array(data.length);
        for (let i = 0, n = data.length; i < n; i++) {
            indices[i] = i;
        }
        jobResult.sortedRows = [];
        for (let i = 0, n = groups.length; i < n; i++) {
            jobResult.sortedRows.push(indices);
        }
    }
    if (jobResult.by == null) {
        jobResult.by = jobResult.fields[0];
        for (let i = 0; i < jobResult.fields.length; i++) {
            if (jobResult.fields[i].toLowerCase().indexOf('score') !== -1) {
                jobResult.by = jobResult.fields[i];
                break;
            }
        }
    }
    if (jobResult.sortByGroup == null) {
        jobResult.sortByGroup = jobResult.groups[0];
    }
    if (jobResult.color == null) {
        jobResult.color = jobResult.fields[0];
    }
    if (jobResult.size == null) {
        jobResult.size = jobResult.fields[0];
    }
    if (jobResult.rowSortOrder == null) {
        jobResult.rowSortOrder = [];
    }
    if (jobResult.rows == null) {
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
    if (jobResult.colorScale == null) {
        let domain = getRange(jobResult.color);
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
    if (jobResult.sizeScaleReversed == null) {
        jobResult.sizeScaleReversed = false;
    }
    if (jobResult.sizeScale == null) {
        let domain = getRange(jobResult.size);
        if (!isNaN(jobResult.options.minSize)) {
            domain[0] = jobResult.options.minSize;
        }
        if (!isNaN(jobResult.options.maxSize)) {
            domain[1] = jobResult.options.maxSize;
        }
        jobResult.sizeScale = scaleLinear().domain(domain).range(jobResult.sizeScaleReversed ? [18, 2] : [2, 18]).clamp(true);
    }
}

export function sortAndFilterJobResult(jobResult) {
    const by = jobResult.by;

    const groups = jobResult.groups;
    const ngroups = groups.length;
    jobResult.sortedRows = [];
    // sort
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

export function sortByGroup(jobResult) {
    const groupIndices = jobResult.sortedFilteredRows[jobResult.groups.indexOf(jobResult.sortByGroup)];
    const indexToRank = [];
    for (let i = 0; i < groupIndices.length; i++) {
        indexToRank[groupIndices[i]] = i + 1;
    }
    jobResult.rows.sort((a, b) => {
        const rankA = indexToRank[a] || a;
        const rankB = indexToRank[b] || b;
        return rankA - rankB;
    });
}

export function updateTopNJobResult(jobResult) {
    let indices = new Set();
    for (let groupIndex = 0; groupIndex < jobResult.groups.length; groupIndex++) {
        const groupIndices = jobResult.sortedFilteredRows[groupIndex];
        const ntop = Math.min(jobResult.ntop, groupIndices.length);
        for (let i = 0; i < ntop; i++) {
            indices.add(groupIndices[i]);
        }
    }
    jobResult.rows = Array.from(indices);
    sortByGroup(jobResult);
}

class JobResultsPanel extends React.PureComponent {


    constructor(props) {
        super(props);
        this.state = {showDialog: false};
    }

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

    getJobResult = (id) => {
        const {jobResults, jobResultId} = this.props;
        if (id == null) {
            id = jobResultId;
        }
        return find(jobResults, item => item.id === id);
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
            searchTokens.push({value: feature, type: X_SEARCH_TOKEN});
            this.props.setTab('embedding');
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

        const {jobResultId, jobResults, classes, searchTokens, tab} = this.props;
        const jobResult = this.getJobResult();
        const selectedFeatures = new Set();
        searchTokens.forEach(token => {
            if (token.type === X_SEARCH_TOKEN) {
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
        if (jobResult != null) {
            // const name = jobResult.name;
            // const params = jobResult.params;
            data = jobResult.data;
            color = jobResult.color;
            by = jobResult.by;
            size = jobResult.size;
            groups = jobResult.groups;
            for (let i = 0; i < groups.length; i++) {
                if (groups[i].length > 2) {
                    rotateHeaders = true;
                    break;
                }
            }
            rows = jobResult.rows;
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
                headerWidth = Math.min(headerWidth, 130);
                headerHeight = Math.cos(45) * headerWidth;
            }
            maxSize = Math.max(jobResult.sizeScale.range()[0], jobResult.sizeScale.range()[1]);
        }

        function formatNumber(val) {
            return (val == null || isNaN(val)) ? "NaN" : numberFormat2f(val);
        }

        const showBrowseJobs = (jobResultId == null && jobResults.length > 0) || jobResults.length > 1;
        return <React.Fragment>
            <Box color="text.primary">
                {showBrowseJobs &&
                <Button style={{marginBottom: 8}} size={"small"} onClick={this.onBrowseJobs} variant="outlined"
                        color="primary">Browse All Results</Button>}
                {jobResult && <React.Fragment>
                    <Typography
                        style={{marginBottom: headerHeight + 8}}
                        component={"h3"}
                        color={"textPrimary"}><b>{jobResult.name}</b>&nbsp;
                        <small>{intFormat(rows.length) + ' / ' + intFormat(jobResult.data.length) + ' features'}</small></Typography>
                    <div style={{paddingTop: 6}}>
                        <Table stickyHeader={true} className={classes.table}>
                            <TableHead>
                                <TableRow>
                                    <TableCell style={{backgroundColor: 'unset', textAlign: 'center'}}
                                               key={'__id'}></TableCell>
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
                                {rows.map(row => {
                                    const id = data[row]['index'];
                                    const selected = selectedFeatures.has(id);
                                    return <TableRow
                                        hover
                                        onClick={(event) => this.handleClick(event, row)}
                                        role="checkbox"
                                        tabIndex={-1}
                                        key={row}
                                    ><TableCell className={classes.rowHeader} component="th"
                                                key={'id'}><Checkbox className={classes.checkbox}
                                                                     checked={selected}/>{id}</TableCell>
                                        {jobResult.columns.map(column => {
                                            const group = groups[column];
                                            const colorField = group + ':' + color;
                                            const sizeField = group + ':' + size;
                                            const byField = group + ':' + by;
                                            const colorValue = data[row][colorField];
                                            const sizeValue = data[row][sizeField];
                                            const byValue = data[row][byField];
                                            const diameter = jobResult.sizeScale(sizeValue);
                                            const title = by + ':' + formatNumber(byValue) + ', '
                                                + color + ':' + formatNumber(colorValue) +
                                                ', ' + size + ':' + formatNumber(sizeValue);

                                            return <TableCell className={classes.td} title={title}
                                                              key={group}>
                                                <div className={classes.dot}
                                                     style={{
                                                         marginLeft: (maxSize - diameter) / 2,
                                                         width: diameter,
                                                         height: diameter,
                                                         backgroundColor: jobResult.colorScale(colorValue)
                                                     }}></div>
                                            </TableCell>;
                                        })}
                                    </TableRow>;
                                })}
                            </TableBody>
                        </Table>
                    </div>
                </React.Fragment>
                }
            </Box>
            <Dialog onClose={this.onCloseJobs} aria-labelledby="job-results-title"
                    open={this.state.showDialog || (tab === 'results' && jobResults.length > 1 && jobResult == null)}>
                <DialogTitle id="job-results-title" onClose={this.onCloseJobs}>
                    Results
                </DialogTitle>
                <DialogContent dividers>
                    <List style={{width: 500}} dense disablePadding component="nav"
                          aria-label="results">
                        {jobResults.map(choice => {
                            let text = choice.name;
                            if (choice.title) {
                                text += ' - ' + choice.title;
                            }
                            return <ListItem alignItems="flex-start"
                                             selected={choice.id === jobResultId}
                                             key={choice.id}
                                             button
                                             onClick={(e) => this.onSelectJob(choice.id)}>
                                <ListItemText
                                    primary={text}
                                    style={{
                                        textOverflow: 'ellipsis',
                                        overflow: 'hidden',
                                        whiteSpace: 'nowrap'
                                    }}/>
                            </ListItem>;
                        })}
                    </List>
                </DialogContent>
            </Dialog>
        </React.Fragment>;
    };
}


const mapStateToProps = state => {
        return {
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

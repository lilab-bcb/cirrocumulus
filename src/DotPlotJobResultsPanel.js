import Box from '@mui/material/Box';
import Typography from '@mui/material/Typography';
import {scaleLinear} from 'd3-scale';
import React, {useState} from 'react';
import {connect} from 'react-redux';
import {setSearchTokens} from './actions';
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
import {Tooltip} from '@mui/material';
import IconButton from '@mui/material/IconButton';
import CloudDownloadIcon from '@mui/icons-material/CloudDownload';
import CirroTooltip from './CirroTooltip';


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

function DotPlotJobResultsPanel(props) {
    const {dataset, jobResult, searchTokens, onSearchTokens} = props;
    const [tip, setTip] = useState({html: ''});

    function onMouseMove(event) {
        setTip({html: event.target.dataset.title, clientX: event.clientX, clientY: event.clientY});
    }

    function onMouseOut(event) {
        setTip({html: ''});
    }


    function exportJobResult(event) {
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
    }


    function toggleAll(event, isSelected) {
        event.stopPropagation();
        isSelected = !isSelected;
        const features = new Set();
        for (let i = 0; i < jobResult.rows.length; i++) {
            const feature = jobResult.data[jobResult.rows[i]]['index'];
            features.add(feature);
        }

        let filteredSearchTokens;
        if (!isSelected) {
            filteredSearchTokens = searchTokens.filter(token => !features.has(token.id));
        } else {
            filteredSearchTokens = searchTokens.slice();
            features.forEach(feature => {
                let found = false;
                for (let i = 0; i < searchTokens.length; i++) {
                    if (searchTokens[i].id === feature) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    filteredSearchTokens.push({
                        id: feature,
                        type: dataset.obs.indexOf(feature) !== -1 ? FEATURE_TYPE.OBS : FEATURE_TYPE.X
                    });
                }
            });
        }
        onSearchTokens(filteredSearchTokens);
    }

    function onRowClick(event, row) {
        event.stopPropagation();
        const feature = jobResult.data[row]['index'];
        let index = -1;
        for (let i = 0; i < searchTokens.length; i++) {
            if (searchTokens[i].id === feature) {
                index = i;
                break;
            }
        }
        if (index === -1) {
            searchTokens.push({
                id: feature,
                type: dataset.obs.indexOf(feature) !== -1 ? FEATURE_TYPE.OBS : FEATURE_TYPE.X
            });
        } else {
            searchTokens.splice(index, 1);
        }
        onSearchTokens(searchTokens.slice());
    }

    function sortColumns() {
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
    }


    const selectedFeatures = new Set();
    searchTokens.forEach(token => {
        if (token.type === FEATURE_TYPE.X) {
            selectedFeatures.add(token.id);
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

    let columns = jobResult.columns;
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


    columns = jobResult.columns.map(column => groups[column]);
    const isRowSelected = (row) => selectedFeatures.has(data[row].index);
    const getRowId = (row) => data[row].index;
    const getTooltip = (row, column) => {
        let title = [];
        tooltipFields.forEach(field => {
            const value = data[row][column + ':' + field];
            title.push(field + ': ' + value);
        });
        title = title.join('<br />');
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

    return <Box color="text.primary">
        <Typography
            style={{marginBottom: headerHeight + 8}}
            component={"h3"}
            color={"textPrimary"}><b>{jobResult.name}</b>&nbsp;
            <small>{intFormat(rows.length) + ' / ' + intFormat(jobResult.data.length) + ' features'}</small></Typography>
        <Tooltip title={"Export"}><IconButton edge={false} size={'small'} aria-label="Export"
                                              onClick={exportJobResult}><CloudDownloadIcon/></IconButton></Tooltip>
        <div style={{paddingTop: 6, position: 'relative'}}>
            {rows.length > 0 &&
                <div style={{position: 'absolute'}}><DotPlotTable isRowSelected={isRowSelected}
                                                                  rows={rows}
                                                                  columns={columns}
                                                                  onMouseMove={onMouseMove}
                                                                  onMouseOut={onMouseOut}
                                                                  getTooltip={getTooltip}
                                                                  sizeScale={jobResult.sizeScale}
                                                                  valueScale={valueScale}
                                                                  rowStart={rowStart}
                                                                  getColor={getColor}
                                                                  getSize={getSize}
                                                                  getRowId={getRowId}
                                                                  columnStart={columnStart}
                                                                  onRowClick={onRowClick}
                                                                  colorScale={jobResult.colorScale}
                                                                  toggleAll={toggleAll}
                                                                  rotateHeaders={rotateHeaders}
                                                                  headerHeight={headerHeight}
                                                                  headerWidth={headerWidth}/>
                    <CirroTooltip html={tip.html} clientX={tip.clientX} clientY={tip.clientY}/></div>}
        </div>
    </Box>;
}


const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        searchTokens: state.searchTokens
    };
};
const mapDispatchToProps = (dispatch) => {
        return {
            onSearchTokens: (payload) => {
                dispatch(setSearchTokens(payload));
            }
        };
    }
;


export default (connect(
    mapStateToProps, mapDispatchToProps
)(DotPlotJobResultsPanel));

import Typography from '@mui/material/Typography';
import {scaleLinear} from 'd3-scale';
import React, {useCallback, useEffect, useState} from 'react';
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
const DotPlotTableMemo = React.memo(DotPlotTable);

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

        if (jobResult.interpolator.scale !== INTERPOLATOR_SCALING_NONE) {
            const domains = [];
            // can scale colors globally, by row, or by column
            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                jobResult.rows.forEach(row => {
                    let min = Number.MAX_VALUE;
                    let max = -Number.MAX_VALUE;
                    jobResult.columns.forEach(column => {
                        const group = jobResult.groups[column];
                        const colorField = group + ':' + jobResult.color;
                        const colorValue = jobResult.data[row][colorField];
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
                        const group = jobResult.groups[column];
                        const colorField = group + ':' + jobResult.color;
                        const colorValue = jobResult.data[row][colorField];
                        if (colorValue != null && !isNaN(colorValue)) {
                            min = colorValue < min ? colorValue : min;
                            max = colorValue > max ? colorValue : max;
                        }
                    });
                    domains.push([min, max]);
                });
            }
            jobResult.valueScale = scaleLinear().range([0, 1]);
            jobResult.domains = domains;
        }
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
    const {dataset, darkMode, jobResult, searchTokens, onSearchTokens} = props;
    const [tip, setTip] = useState({html: ''});
    const [headerSize, setHeaderSize] = useState({width: 0, height: 0});
    const [rotateHeaders, setRotateHeaders] = useState(false);
    const [selectedFeatures, setSelectedFeatures] = useState(new Set());
    const [tooltipFields, setTooltipFields] = useState([]);

    const [columns, setColumns] = useState([]);

    useEffect(() => {
        let _rotateHeaders = false;
        for (let i = 0; i < jobResult.groups.length; i++) {
            if (jobResult.groups[i].length > 2) {
                _rotateHeaders = true;
                break;
            }
        }
        let headerWidth = 0;
        let headerHeight = 12;
        if (_rotateHeaders) {

            const d = document.createElement('span');
            d.style.position = 'absolute';
            d.style.left = '-1000000px';
            d.className = '.MuiTableCell-head .MuiTableCell-root';
            document.body.append(d);
            for (let i = 0; i < jobResult.groups.length; i++) {
                d.innerText = jobResult.groups[i];
                headerWidth = Math.max(headerWidth, 2 + d.getBoundingClientRect().width);
            }
            d.remove();
            headerWidth += 2; // prevent overflow
            headerWidth = Math.min(headerWidth, 300);
            headerHeight = Math.cos(45) * headerWidth;
            setHeaderSize({width: headerWidth, height: headerHeight});
        }
        setHeaderSize({width: headerWidth, height: headerHeight});
        setRotateHeaders(_rotateHeaders);
    }, [jobResult.groups]);

    const rowStart = useCallback(
        (row, rowIndex) => {
            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_FEATURE) {
                jobResult.valueScale.domain(jobResult.domains[rowIndex]);
            }
        },
        [jobResult.interpolator, jobResult.valueScale, jobResult.domains]
    );
    useEffect(() => {
        const selection = new Set();
        searchTokens.forEach(token => {
            if (token.type === FEATURE_TYPE.X) {
                selection.add(token.id);
            }
        });
        setSelectedFeatures(selection);

    }, [searchTokens]);
    useEffect(() => {
        const fields = [jobResult.by];
        if (jobResult.color !== jobResult.by) {
            fields.push(jobResult.color);
        }
        const isSizeScaled = jobResult.size !== 'none';
        if (jobResult.size !== jobResult.by && isSizeScaled && jobResult.size !== jobResult.color) {
            fields.push(jobResult.size);
        }

        jobResult.fields.forEach(field => {
            if (fields.indexOf(field) === -1) {
                fields.push(field);
            }
        });
        setTooltipFields(fields);
    }, [jobResult.color, jobResult.by, jobResult.size, jobResult.fields]);


    useEffect(() => {
        setColumns(jobResult.columns.map(c => jobResult.groups[c]));
    }, [jobResult.columns, jobResult.groups]);
    const isRowSelected = useCallback(
        (row) => {
            return selectedFeatures.has(jobResult.data[row].index);
        },
        [selectedFeatures, jobResult.data]
    );

    const getRowId = useCallback(
        (row) => {
            return jobResult.data[row].index;
        },
        [jobResult.data]
    );

    const getTooltip = useCallback(
        (row, column) => {
            let title = [];
            const item = jobResult.data[row];
            tooltipFields.forEach(field => {
                const key = column + ':' + field;
                const value = item[key];
                if (value != null) {
                    title.push(field + ': ' + value);
                }
            });
            title = title.join('<br />');
            return title;
        },
        [jobResult.data, tooltipFields]
    );


    const getColor = useCallback(
        (row, column) => {
            return jobResult.data[row][column + ':' + jobResult.color];
        },
        [jobResult.data, jobResult.color]
    );

    const getSize = useCallback(
        (row, column) => {
            return jobResult.data[row][column + ':' + jobResult.size];
        },
        [jobResult.data, jobResult.size]
    );


    const columnStart = useCallback(
        (column, columnIndex) => {
            if (jobResult.interpolator.scale === INTERPOLATOR_SCALING_MIN_MAX_CATEGORY) {
                jobResult.valueScale.domain(jobResult.domains[columnIndex]);
            }
        },
        [jobResult.interpolator, jobResult.valueScale, jobResult.domains]
    );

    const onMouseMove = useCallback(
        (event) => {
            setTip({
                html: event.target.dataset.title == null ? '' : event.target.dataset.title,
                clientX: event.clientX,
                clientY: event.clientY
            });
        }, []
    );

    const onMouseOut = useCallback(
        (event) => {
            setTip({html: ''});
        }, []
    );

    const toggleAll = useCallback(
        (event, isSelected) => {
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
        }, [jobResult.rows, jobResult.data, searchTokens]
    );

    const onRowClick = useCallback(
        (event, row) => {
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
        }, [jobResult.data, searchTokens]
    );


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


    // function sortColumns() {
    //     const sortOrder = jobResult.columnSortOrder;
    //     const groups = jobResult.groups;
    //     const data = jobResult.data;
    //     const indices = new Array(groups.length);
    //     for (let i = 0, n = groups.length; i < n; i++) {
    //         indices[i] = i;
    //     }
    //     const nsortFields = sortOrder.length;
    //     if (nsortFields.length > 0) { // sort by gene
    //         const sortFieldRows = [];
    //         for (let i = 0; i < nsortFields; i++) {
    //             const field = sortOrder[i].field;
    //             let index;
    //             for (let j = 0; j < data.length; j++) {
    //                 if (data[j]['index'] === field) {
    //                     index = j;
    //                     break;
    //                 }
    //             }
    //             sortFieldRows.push(index);
    //         }
    //
    //         indices.sort((a, b) => {
    //             for (let i = 0; i < nsortFields; i++) {
    //                 const row = sortFieldRows[i];
    //                 const val = data[row][jobResult.color] - data[row][jobResult.color];
    //                 if (val !== 0) {
    //                     return val;
    //                 }
    //                 const val2 = data[row][jobResult.size] - data[row][jobResult.size];
    //                 if (val2 !== 0) {
    //                     return val2;
    //                 }
    //             }
    //             return a - b;
    //         });
    //     }
    //     jobResult.columns = indices;
    // }


    return <div>
        <Typography
            style={{marginBottom: headerSize.height + 8}}
            component={"h3"}
            color={"textPrimary"}><b>{jobResult.name}</b>&nbsp;
            <small>{intFormat(jobResult.rows.length) + ' / ' + intFormat(jobResult.data.length) + ' features'}</small></Typography>
        <Tooltip title={"Export"}><IconButton edge={false} size={'small'} aria-label="Export"
                                              onClick={exportJobResult}><CloudDownloadIcon/></IconButton></Tooltip>
        <div style={{paddingTop: 6, position: 'relative'}}>
            {jobResult.rows.length > 0 &&
                <div style={{position: 'absolute', width: '100%'}}>
                    <DotPlotTableMemo isRowSelected={isRowSelected}
                                      rows={jobResult.rows}
                                      columns={columns}
                                      colorScale={jobResult.colorScale}
                                      darkMode={darkMode}
                                      rotateHeaders={rotateHeaders}
                                      headerHeight={headerSize.height}
                                      headerWidth={headerSize.width}
                                      onMouseMove={onMouseMove}
                                      onMouseOut={onMouseOut}
                                      getTooltip={getTooltip}
                                      sizeScale={jobResult.sizeScale}
                                      valueScale={jobResult.valueScale}
                                      rowStart={rowStart}
                                      getColor={getColor}
                                      getSize={getSize}
                                      getRowId={getRowId}
                                      columnStart={columnStart}
                                      onRowClick={onRowClick}
                                      toggleAll={toggleAll}/>
                    <CirroTooltip html={tip.html} clientX={tip.clientX} clientY={tip.clientY}/></div>}
        </div>
    </div>;
}


const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        searchTokens: state.searchTokens,
        darkMode: state.chartOptions.darkMode
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

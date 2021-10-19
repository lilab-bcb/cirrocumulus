import {Checkbox, Table, TableBody, TableCell, TableHead, TableRow} from '@mui/material';
import Box from '@mui/material/Box';
import withStyles from '@mui/styles/withStyles';
import React from 'react';

const styles = theme => ({
    dot: {
        borderRadius: '50%',
        display: 'inlineBlock',
        border: '1px solid lightgray'
    },
    table: {
        borderCollapse: 'collapse',
        position: 'relative',
        width: 'unset'
    },
    tr: {
        cursor: 'pointer'
    },
    deleteTr: {
        cursor: 'pointer',
        '& span': {
            display: 'none'
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
        width: '100%'
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
        whiteSpace: 'nowrap'
    }
});


function DotPlotTable(props) {
    const {
        classes,
        onRowClick,
        onMouseMove,
        onMouseOut,
        toggleAll,
        rows,
        columns,
        isRowSelected,
        rowStart,
        columnStart,
        getRowId,
        getColor,
        getSize,
        getTooltip,
        sizeScale,
        valueScale,
        colorScale,
        headerWidth,
        rotateHeaders
    } = props;

    let selectAllChecked = true;
    const maxSize = Math.max(sizeScale.range()[0], sizeScale.range()[1]);
    for (let i = 0; i < rows.length; i++) {
        const selected = isRowSelected(rows[i]);
        if (!selected) {
            selectAllChecked = false;
            break;
        }
    }

    return <div data-testid={'dot-plot-table'}>
        <Box color="text.primary">
            <div>
                {rows.length > 0 &&
                <Table onMouseMove={onMouseMove} onMouseOut={onMouseOut} stickyHeader={true}
                       className={classes.table}>
                    <TableHead>
                        <TableRow>
                            <TableCell className={classes.rowHeader} component="th"
                                       style={{backgroundColor: 'unset', textAlign: 'left'}}
                                       key={'__id'}><Checkbox
                                onClick={(event) => toggleAll(event, selectAllChecked)}
                                className={classes.checkbox}
                                checked={selectAllChecked}/></TableCell>
                            {columns.map(column => {
                                if (rotateHeaders) {
                                    return <TableCell
                                        title={column}
                                        className={classes.rotateHeader}
                                        key={column}>
                                        <div className={classes.rotateHeaderDiv}><span
                                            style={{width: headerWidth}}
                                            className={classes.rotateHeaderSpan}>{column}</span></div>
                                    </TableCell>;
                                } else {
                                    return <TableCell style={{backgroundColor: 'unset'}}
                                                      className={classes.td}
                                                      key={column}>{column}</TableCell>;
                                }
                            })}
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {rows.map((row, rowIndex) => {
                            rowStart(row, rowIndex);
                            const id = getRowId(row);
                            const selected = isRowSelected(row);
                            return <TableRow
                                className={classes.tr}
                                hover
                                onClick={(event) => onRowClick(event, row)}
                                role="checkbox"
                                tabIndex={-1}
                                key={row}
                            ><TableCell className={classes.rowHeader} component="th"
                                        key={'id'}><Checkbox className={classes.checkbox}
                                                             checked={selected}/>{id}</TableCell>
                                {columns.map((column, columnIndex) => {
                                    columnStart(column, columnIndex);
                                    const colorValue = getColor(row, column);
                                    const sizeValue = getSize(row, column);
                                    const title = getTooltip(row, column);
                                    if (colorValue == null || isNaN(colorValue)) {
                                        return <TableCell className={classes.td} data-title={title}
                                                          key={column}/>;
                                    }
                                    const diameter = sizeScale(sizeValue);
                                    const colorValueScaled = valueScale ? valueScale(colorValue) : colorValue;
                                    const backgroundColor = colorScale(colorValueScaled);
                                    return <TableCell className={classes.td} data-title={title}
                                                      key={column}>
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
        </Box>
    </div>;
}

export default withStyles(styles)((DotPlotTable));


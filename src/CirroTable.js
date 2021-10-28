import {
    Checkbox,
    ListItemText,
    Select,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TablePagination,
    TableRow,
    TableSortLabel
} from '@mui/material';
import Box from '@mui/material/Box';
import {visuallyHidden} from '@mui/utils';
import IconButton from '@mui/material/IconButton';
import {intFormat} from './formatters';
import React, {useState} from 'react';
import {NATSORT} from './util';
import TextField from '@mui/material/TextField';
import InputAdornment from '@mui/material/InputAdornment';
import ClearIcon from '@mui/icons-material/Clear';
import FormControl from '@mui/material/FormControl';
import MenuItem from '@mui/material/MenuItem';
import {findIndex} from 'lodash';

export default function CirroTable(props) {
    const {isSelected, onItemClick, rows, columns, rowId, renderCell, onSearchText, searchText} = props;

    const [rowsPerPage, setRowsPerPage] = useState(30);
    const [page, setPage] = useState(0);
    const [order, setOrder] = useState('asc');
    const [orderBy, setOrderBy] = useState(null);

    let filteredRows = rows;
    let searchTextLower = searchText.toLowerCase().trim();
    const allColumns = columns;
    const visibleColumns = allColumns.filter(column => column.visible !== false);
    const [forceUpdate, setForceUpdate] = useState(false);

    if (searchTextLower != '') {
        let searchColumns = visibleColumns;
        const sepIndex = searchTextLower.indexOf(':'); // field search
        if (sepIndex !== -1) {
            let field = searchTextLower.substring(0, sepIndex);
            for (let i = 0; i < allColumns.length; i++) {
                if (field == allColumns[i].label.toLowerCase().replace(' ', '_')) {
                    searchColumns = [allColumns[i]];
                    searchTextLower = searchTextLower.substring(sepIndex + 1).trim();
                    break;
                }
            }
        }
        if (searchTextLower != '') { // need to check again as searchTextLower might have been updated
            const ncolumns = searchColumns.length;
            const matches = (value) => value != null && ('' + value).toLowerCase().indexOf(searchTextLower) !== -1;
            filteredRows = filteredRows.filter(choice => {
                for (let i = 0; i < ncolumns; i++) {
                    const value = choice[searchColumns[i].field];
                    if (matches(value)) {
                        return true;
                    }
                }
                return false;
            });
        }
    }
    if (orderBy != null) {
        filteredRows.sort((item1, item2) => {
            return NATSORT(item1[orderBy], item2[orderBy]);
        });
        if (order == 'desc') {
            filteredRows.reverse();
        }
    }
    const handleRequestSort = (event, property) => {
        const isAsc = orderBy === property && order === 'asc';
        setOrder(isAsc ? 'desc' : 'asc');
        setOrderBy(property);
    };

    const handleChangePage = (event, newPage) => {
        setPage(newPage);
    };


    function handleItemClick(item) {
        onItemClick(item);

    }

    function handleColumnsChange(event) {
        const value = event.target.value;
        columns.forEach(column => {
            column.visible = value.indexOf(column.field) !== -1;
        });
        setForceUpdate(!forceUpdate); // TODO redux style
    }

    function handleClearSearchText() {
        setPage(0);
        onSearchText('');
    }

    function onSearchChange(event) {
        setPage(0);
        onSearchText(event.target.value);
    }

    return <div>
        <TextField size="small" style={{paddingTop: 6}} type="text" placeholder={"Search"}
                   value={searchText}
                   sx={{width: '80%', maxWidth: 800}}
                   onChange={onSearchChange}
                   helperText={"Search a specific field by typing the field name followed by a colon \":\" and then the term you are looking for"}
                   InputProps={searchText.trim() !== '' ? {
                       endAdornment:
                           <InputAdornment position="end">
                               <IconButton aria-label="clear" onClick={handleClearSearchText}
                                           size="large">
                                   <ClearIcon/>
                               </IconButton>
                           </InputAdornment>
                   } : null}
        />
        <FormControl sx={{m: 1}}>
            <Select
                size={"small"}
                variant={'standard'}
                labelId="selector-columns-label"
                id="selector-columns"
                multiple
                value={visibleColumns.map(item => item.field)}
                onChange={handleColumnsChange}
                renderValue={(selected) => 'Column Visibility'}
            >
                {allColumns.map((item) => (
                    <MenuItem key={item.field} value={item.field}>
                        <Checkbox checked={findIndex(visibleColumns, c => c.field === item.field) !== -1}/>
                        <ListItemText primary={item.label}/>
                    </MenuItem>
                ))}
            </Select>
        </FormControl>
        <TableContainer>
            <Table stickyHeader size={"small"} padding={"normal"}>
                <TableHead>
                    <TableRow>
                        {visibleColumns.map((column) => (
                            <TableCell
                                key={column.field}
                                align={column.numeric ? 'right' : 'left'}
                                padding={column.disablePadding ? 'none' : 'normal'}
                                sortDirection={orderBy === column.field ? order : false}
                            >
                                <TableSortLabel
                                    active={orderBy === column.field}
                                    direction={orderBy === column.field ? order : 'asc'}
                                    onClick={event => handleRequestSort(event, column.field)}
                                >
                                    {column.label}
                                    {orderBy === column.field ? (
                                        <Box component="span" sx={visuallyHidden}>
                                            {order === 'desc' ? 'sorted descending' : 'sorted ascending'}
                                        </Box>
                                    ) : null}
                                </TableSortLabel>
                            </TableCell>
                        ))}
                    </TableRow>
                </TableHead>
                <TableBody>
                    {filteredRows
                        .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                        .map((row, rowIndex) => {
                            return (
                                <TableRow selected={isSelected(row, rowIndex)} hover role="checkbox"
                                          tabIndex={-1} key={rowId(row, rowIndex)}
                                          onClick={(e) => handleItemClick(row)}>
                                    {visibleColumns.map((column, columnIndex) => {
                                        return (
                                            <TableCell key={column.field} align={column.align}>
                                                {renderCell(row, column, columnIndex)}
                                            </TableCell>
                                        );
                                    })}
                                </TableRow>
                            );
                        })}
                </TableBody>
            </Table>
        </TableContainer>
        <TablePagination
            labelDisplayedRows={({from, to, count}) => {
                from = intFormat(from);
                to = intFormat(to);
                count = intFormat(count);
                return from + '-' + to + ' of ' + count + (rows.length !== filteredRows.length ? ' (filtered from ' + intFormat(rows.length) + ' total)' : '');
            }}
            rowsPerPageOptions={[50]}
            component="div"
            count={filteredRows.length}
            rowsPerPage={rowsPerPage}
            page={page}
            onPageChange={handleChangePage}
        />
    </div>;
}
import {
    Divider,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TablePagination,
    TableRow,
    TableSortLabel,
    Tooltip,
    Typography
} from '@mui/material';
import InfoIcon from '@mui/icons-material/Info';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogContent from '@mui/material/DialogContent';
import IconButton from '@mui/material/IconButton';
import InputAdornment from '@mui/material/InputAdornment';
import Popover from '@mui/material/Popover';
import TextField from '@mui/material/TextField';
import ClearIcon from '@mui/icons-material/Clear';
import FolderOpenIcon from '@mui/icons-material/FolderOpen';
import ReactMarkdown from 'markdown-to-jsx';
import React, {useState} from 'react';
import {connect} from 'react-redux';
import {OPEN_DATASET_DIALOG, setDialog} from './actions';
import {NATSORT, REACT_MD_OVERRIDES} from './util';
import {intFormat} from './formatters';
import {visuallyHidden} from '@mui/utils';
import Box from '@mui/material/Box';

import {find} from 'lodash';

export function DatasetSelector(props) {
    const [searchText, setSearchText] = useState('');
    const [datasetDetailsEl, setDatasetDetailsEl] = useState(null);
    const [selectedDataset, setSelectedDataset] = useState(null);
    const [rowsPerPage, setRowsPerPage] = useState(30);
    const [page, setPage] = useState(0);
    const [order, setOrder] = useState('asc');
    const [orderBy, setOrderBy] = useState(null);
    const {dataset, datasetChoices, dialog} = props;
    const handleRequestSort = (event, property) => {
        const isAsc = orderBy === property && order === 'asc';
        setOrder(isAsc ? 'desc' : 'asc');
        setOrderBy(property);
    };

    const handleChangePage = (event, newPage) => {
        setPage(newPage);
    };

    function handleClick(event) {
        props.handleDialog(OPEN_DATASET_DIALOG);
        setSearchText('');
    }

    function handleClose() {
        props.handleDialog(null);
        setSearchText('');
    }

    function handleCloseDatasetDetails(event) {
        setDatasetDetailsEl(null);
        setSelectedDataset(null);
    }

    function handleListItemDetailsClick(event, id) {
        event.stopPropagation();
        setDatasetDetailsEl(event.currentTarget);
        setSelectedDataset(find(datasetChoices, item => id === item.id));
    }

    function handleListItemClick(id) {
        const selectedId = props.dataset != null ? props.dataset.id : null;
        if (id !== selectedId) {
            props.onChange(id);
        }
        props.handleDialog(null);
        setSearchText('');
    }

    function handleClearSearchText() {
        setPage(0);
        setSearchText('');
    }

    function onSearchChange(event) {
        setPage(0);
        setSearchText(event.target.value);
    }


    const selectedId = dataset != null ? dataset.id : null;

    if (datasetChoices.length <= 1 && selectedId != null) {
        return null;
    }
    const open = dialog === OPEN_DATASET_DIALOG;

    let filteredChoices = datasetChoices;
    let searchTextLower = searchText.toLowerCase().trim();
    const allColumns = props.serverInfo.datasetSelectorColumns || [{id: 'name', label: 'Name'}, {
        id: 'species',
        label: 'Species'
    }, {
        id: 'title',
        label: 'Title'
    }];
    const visibleColumns = allColumns.filter(column => column.visible !== false);

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
            filteredChoices = filteredChoices.filter(choice => {
                for (let i = 0; i < ncolumns; i++) {
                    const value = choice[searchColumns[i].id];
                    if (value != null && ('' + value).toLowerCase().indexOf(searchTextLower) !== -1) {
                        return true;
                    }
                }
                return false;
            });
        }
    }
    if (orderBy != null) {
        filteredChoices.sort((item1, item2) => {
            return NATSORT(item1[orderBy], item2[orderBy]);
        });
        if (order == 'desc') {
            filteredChoices.reverse();
        }
    }
    const datasetDetailsOpen = Boolean(datasetDetailsEl);
    // const group2Items = groupBy(filteredChoices, item => item[datasetSelectorGroupBy] || '');
    // const groups = Object.keys(group2Items);
    // groups.sort(NATSORT);
    // const otherIndex = groups.indexOf("");
    // if (otherIndex !== -1) { // move Other (empty string) to end
    //     groups.push(groups.splice(otherIndex, 1)[0]);
    // }

    return (
        <React.Fragment>
            <Popover
                id={"dataset-details-selector"}
                open={datasetDetailsOpen}
                anchorEl={datasetDetailsEl}
                onClose={handleCloseDatasetDetails}
                anchorOrigin={{
                    vertical: 'bottom',
                    horizontal: 'center'
                }}
                transformOrigin={{
                    vertical: 'top',
                    horizontal: 'center'
                }}
            >
                <div style={{width: 500}}>
                    {selectedDataset && selectedDataset.description &&
                    <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                   children={selectedDataset.description}/>}
                </div>
            </Popover>
            {selectedId == null && <Button variant="contained" onClick={handleClick}
                                           color="primary" startIcon={<FolderOpenIcon/>}>Open</Button>}
            {selectedId != null &&
            <Tooltip title={'Open'}><IconButton onClick={handleClick}
                                                size="large"><FolderOpenIcon/></IconButton></Tooltip>}

            <Dialog
                fullWidth={true}
                open={open}
                onClose={handleClose}
                maxWidth={"xl"}
            >
                <DialogContent sx={{height: '100vh'}}>
                    <TextField size="small" style={{paddingTop: 6}} type="text" placeholder={"Search Datasets"}
                               value={searchText}
                               onChange={onSearchChange}
                               fullWidth={true}
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
                    <div>
                        <Typography
                            variant="subtitle1">{(filteredChoices.length != datasetChoices.length ? (intFormat(filteredChoices.length) + ' / ') : '') + intFormat(datasetChoices.length) + ' Datasets'}</Typography>
                        <Divider/>
                        <TableContainer>
                            <Table stickyHeader size={"small"} padding={"normal"}>
                                <TableHead>
                                    <TableRow>
                                        {visibleColumns.map((column) => (
                                            <TableCell
                                                key={column.id}
                                                align={column.numeric ? 'right' : 'left'}
                                                padding={column.disablePadding ? 'none' : 'normal'}
                                                sortDirection={orderBy === column.id ? order : false}
                                            >
                                                <TableSortLabel
                                                    active={orderBy === column.id}
                                                    direction={orderBy === column.id ? order : 'asc'}
                                                    onClick={event => handleRequestSort(event, column.id)}
                                                >
                                                    {column.label}
                                                    {orderBy === column.id ? (
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
                                    {filteredChoices
                                        .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                                        .map((row) => {
                                            return (
                                                <TableRow selected={row.id === selectedId} hover role="checkbox"
                                                          tabIndex={-1} key={row.id}
                                                          onClick={(e) => handleListItemClick(row.id)}>
                                                    {visibleColumns.map((column, columnIndex) => {
                                                        const value = row[column.id];
                                                        return (
                                                            <TableCell key={column.id} align={column.align}>
                                                                {column.format && typeof value === 'number'
                                                                    ? column.format(value)
                                                                    : value}
                                                                {columnIndex === 0 && row.description != null && row.description !== '' &&
                                                                <IconButton
                                                                    onClick={(e) => handleListItemDetailsClick(e, row.id)}
                                                                    edge="end"
                                                                    aria-label="description"
                                                                    size="small">
                                                                    <InfoIcon/>
                                                                </IconButton>}
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
                            rowsPerPageOptions={[30]}
                            component="div"
                            count={filteredChoices.length}
                            rowsPerPage={rowsPerPage}
                            page={page}
                            onPageChange={handleChangePage}

                        />
                        {/*{groups.map(species => {*/}
                        {/*    const speciesText = species === '' ? 'Other' : species;*/}
                        {/*    const choices = group2Items[species];*/}
                        {/*    choices.sort((a, b) => NATSORT(a.name, b.name));*/}
                        {/*    return (*/}
                        {/*        <React.Fragment key={species}>*/}
                        {/*            <Typography component={"h2"}>{speciesText}</Typography>*/}
                        {/*            <List dense disablePadding component="nav">*/}
                        {/*                {choices.map(choice => {*/}
                        {/*                    let text = choice.name;*/}
                        {/*                    if (choice.title) {*/}
                        {/*                        text += ' - ' + choice.title;*/}
                        {/*                    }*/}
                        {/*                    return (*/}
                        {/*                        <ListItem alignItems="flex-start"*/}
                        {/*                                  selected={choice.id === selectedId}*/}
                        {/*                                  key={choice.id}*/}
                        {/*                                  button*/}
                        {/*                                  onClick={(e) => handleListItemClick(choice.id)}>*/}
                        {/*                            <ListItemText*/}
                        {/*                                primary={text}*/}
                        {/*                                style={{*/}
                        {/*                                    textOverflow: 'ellipsis',*/}
                        {/*                                    overflow: 'hidden',*/}
                        {/*                                    whiteSpace: 'nowrap'*/}
                        {/*                                }}/>*/}
                        {/*                            <ListItemSecondaryAction>*/}
                        {/*                                <IconButton*/}
                        {/*                                    onClick={(e) => handleListItemDetailsClick(e, choice.id)}*/}
                        {/*                                    edge="end"*/}
                        {/*                                    aria-label="summary"*/}
                        {/*                                    size="large">*/}
                        {/*                                    <InfoIcon/>*/}
                        {/*                                </IconButton>*/}
                        {/*                            </ListItemSecondaryAction>*/}
                        {/*                        </ListItem>*/}
                        {/*                    );*/}
                        {/*                })}*/}
                        {/*            </List></React.Fragment>*/}
                        {/*    );*/}
                        {/*})}*/}
                    </div>
                </DialogContent>
            </Dialog>
        </React.Fragment>
    );

}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        dialog: state.dialog,
        datasetChoices: state.datasetChoices,
        serverInfo: state.serverInfo
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleDialog: value => {
            dispatch(setDialog(value));
        }

    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(DatasetSelector));


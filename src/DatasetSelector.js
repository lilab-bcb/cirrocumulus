import {Divider, Tooltip} from '@mui/material';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogContent from '@mui/material/DialogContent';
import IconButton from '@mui/material/IconButton';
import Popover from '@mui/material/Popover';
import FolderOpenIcon from '@mui/icons-material/FolderOpen';
import ReactMarkdown from 'markdown-to-jsx';
import React, {useState} from 'react';
import {connect} from 'react-redux';
import {OPEN_DATASET_DIALOG, setDialog} from './actions';
import {find} from 'lodash';
import {REACT_MD_OVERRIDES} from './util';
import CirroTable from './CirroTable';
import InfoIcon from '@mui/icons-material/Info';
import Box from '@mui/material/Box';
import Typography from '@mui/material/Typography';
import {DATASET_FIELDS} from './EditNewDatasetDialog';

export function DatasetSelector(props) {

    const [datasetDetailsEl, setDatasetDetailsEl] = useState(null);
    const [selectedDataset, setSelectedDataset] = useState(null);
    const [searchText, setSearchText] = useState('');
    const {dataset, datasetChoices, datasetSelectorColumns, dialog, handleDialog, onChange} = props;


    function handleCloseDatasetDetails(event) {
        setDatasetDetailsEl(null);
        setSelectedDataset(null);
    }

    function onSearchText(text) {
        setSearchText(text);
    }

    function renderCell(item, column, columnIndex) {
        let value = item[column.field];
        value = column.format && typeof value === 'number'
            ? column.format(value)
            : value;
        return <>{value}
            {columnIndex === 0 && item.description != null && item.description !== '' && <IconButton
                onClick={(e) => handleListItemDetailsClick(e, item.id)}
                edge="end"
                aria-label="description"
                size="small">
                <InfoIcon/>
            </IconButton>}
        </>;
    }

    function handleListItemDetailsClick(event, id) {
        event.stopPropagation();
        setDatasetDetailsEl(event.currentTarget);
        setSelectedDataset(find(datasetChoices, item => id === item.id));
    }

    function isSelected(item) {
        const selectedId = dataset != null ? dataset.id : null;
        return item.id === selectedId;
    }

    function onItemClick(item) {
        const selectedId = dataset != null ? dataset.id : null;
        const id = item.id;
        if (id !== selectedId) {
            onChange(id);
        }
        handleDialog(null);
        setSearchText('');
    }


    function handleClick(event) {
        handleDialog(OPEN_DATASET_DIALOG);
        setSearchText('');
    }

    function handleClose() {
        handleDialog(null);
        setSearchText('');
    }

    const selectedId = dataset != null ? dataset.id : null;

    if (datasetChoices.length <= 1 && selectedId != null) {
        return null;
    }
    const open = dialog === OPEN_DATASET_DIALOG;
    const datasetDetailsOpen = Boolean(datasetDetailsEl);
    return (
        <>
            {selectedDataset && <Popover
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
                <Box style={{width: 500, padding: '1em'}}>
                    <Typography variant="h6">{selectedDataset.name}</Typography>
                    {DATASET_FIELDS.filter(item => selectedDataset[item.fieldName]).map(item => <div
                        key={item.fieldName}><Divider/>
                        <Typography variant={"subtitle2"}>{item.label}</Typography>{item.fieldName !== 'description' &&
                        <Typography variant="body2"> {selectedDataset[item.fieldName]}</Typography>}
                        {item.fieldName === 'description' &&
                        <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                       children={selectedDataset[item.fieldName]}/>}</div>)}
                </Box>
            </Popover>}
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
                    <CirroTable rows={datasetChoices} columns={datasetSelectorColumns} onItemClick={onItemClick}
                                isSelected={isSelected} renderCell={renderCell} onSearchText={onSearchText}
                                searchText={searchText} rowId={item => item.id}/>
                </DialogContent>
            </Dialog>
        </>
    );

}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        dialog: state.dialog,
        datasetChoices: state.datasetChoices,
        datasetSelectorColumns: state.serverInfo.datasetSelectorColumns
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


import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';

import DialogTitle from '@mui/material/DialogTitle';
import React, {useState} from 'react';
import {connect} from 'react-redux';
import {deleteDataset, setDialog} from './actions';

function DeleteDatasetDialog(props) {
    const {dataset, handleCancel, handleDelete} = props;
    const [loading, setLoading] = useState(false);

    function onClose() {
        handleCancel();
    }

    function onDelete() {
        setLoading(true);
        handleDelete({dataset: dataset});
    }


    return (
        <Dialog
            open={true}
            onClose={onClose}
            aria-labelledby="delete-dataset-dialog-title"
            fullWidth={true}
            maxWidth={'sm'}
        >
            <DialogTitle id="delete-dataset-dialog-title">Delete Dataset</DialogTitle>
            <DialogContent>
                {dataset && <h3>Are you sure you want to delete {dataset.name}?</h3>}
            </DialogContent>
            <DialogActions>
                <Button disabled={loading} onClick={onClose}>
                    Cancel
                </Button>
                <Button disabled={loading} onClick={onDelete} variant="contained"
                        color="primary">
                    Delete
                </Button>
            </DialogActions>
        </Dialog>
    );

}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        email: state.email
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleDelete: value => {
            dispatch(deleteDataset(value));
        },

        handleCancel: value => {
            dispatch(setDialog(null));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(DeleteDatasetDialog));


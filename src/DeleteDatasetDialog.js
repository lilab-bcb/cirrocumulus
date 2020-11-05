import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import React from 'react';
import {connect} from 'react-redux';
import {deleteDataset, setDialog} from './actions';

class DeleteDatasetDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            open: true,
        };
    }

    handleClose = () => {
        this.props.handleCancel();
    };

    handleDelete = () => {
        this.setState({loading: true});
        this.props.handleDelete({dataset: this.props.dataset});
    };

    render() {
        return (
            <Dialog
                open={this.state.open}
                onClose={this.handleClose}
                aria-labelledby="delete-dataset-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="delete-dataset-dialog-title">Delete Dataset</DialogTitle>
                <DialogContent>
                    {this.props.dataset && <h3>Are you sure you want to delete {this.props.dataset.name}?</h3>}
                </DialogContent>
                <DialogActions>
                    <Button disabled={this.state.loading} onClick={this.handleClose} color="primary">
                        Cancel
                    </Button>
                    <Button disabled={this.state.loading} onClick={this.handleDelete} color="primary">
                        Delete
                    </Button>
                </DialogActions>
            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        email: state.email,
    };
};
const mapDispatchToProps = dispatch => {
    return {

        handleDelete: value => {
            dispatch(deleteDataset(value));
        },

        handleCancel: value => {
            dispatch(setDialog(null));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DeleteDatasetDialog));


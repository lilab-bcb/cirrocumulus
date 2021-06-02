import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {saveView, setDialog, setMessage} from './actions';

class SaveDatasetViewDialog extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {
            name: ''
        };
    }

    onNameChanged = (event) => {
        this.setState({name: event.target.value});
    };

    handleClose = () => {
        this.props.handleCancel();
    };

    handleSave = () => {
        let name = this.state.name.trim();
        this.props.handleSave({name: name});
    };


    render() {
        const {name} = this.state;
        const create = true;
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="edit-dataset-view-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="edit-dataset-view-dialog-title">{create
                    ? 'Save'
                    : 'Edit'} Link</DialogTitle>
                <DialogContent>
                    <TextField
                        required={true}
                        value={name}
                        autoComplete="off"
                        onChange={this.onNameChanged}
                        margin="dense"
                        label="Name"
                        fullWidth
                    />

                </DialogContent>
                <DialogActions>
                    <Button onClick={this.handleClose}>
                        Cancel
                    </Button>
                    <Button disabled={name.trim().length === 0} onClick={this.handleSave} variant="contained"
                            color="primary">
                        Save
                    </Button>
                </DialogActions>
            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {};
};
const mapDispatchToProps = dispatch => {
    return {
        handleSave: value => {
            dispatch(saveView(value));
        },
        handleCancel: value => {
            dispatch(setDialog(null));
        },
        handleError: value => {
            dispatch(setMessage(value));
        },

    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(SaveDatasetViewDialog));


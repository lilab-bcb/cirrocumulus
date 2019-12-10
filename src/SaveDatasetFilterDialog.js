import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {saveDatasetFilter, setDialog, setMessage} from './actions';

class SaveDatasetFilterDialog extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {
            open: true,
            name: '',
            notes: ''
        };
    }

    componentDidMount() {
        if (this.props.savedFilter == null) {
            this.setState({name: '', notes: '', create: true});
        } else {
            this.setState({name: this.props.savedFilter.name, notes: this.props.savedFilter.notes});
        }
    }

    onNameChanged = (event) => {
        this.setState({name: event.target.value});
    };

    onNotesChanged = (event) => {
        this.setState({notes: event.target.value});
    };

    handleClose = () => {
        this.props.handleCancel();
    };

    handleSave = () => {
        let name = this.state.name.trim();
        let notes = this.state.notes;
        if (notes != null) {
            notes = notes.trim();
        }
        this.setState({loading: true});
        this.props.handleSave({name: name, notes: notes});
    };


    render() {
        let {name, notes, create} = this.state;
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="edit-dataset-filter-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="edit-dataset-filter-dialog-title">{create
                    ? 'Create'
                    : 'Edit'} Filter</DialogTitle>
                <DialogContent>

                    <TextField
                        required={true}
                        value={name}
                        onChange={this.onNameChanged}
                        margin="dense"
                        label="Name"
                        fullWidth
                    />

                    <TextField
                        required={false}
                        value={notes}
                        onChange={this.onNotesChanged}
                        margin="dense"
                        label="Notes"
                        fullWidth
                    />

                </DialogContent>
                <DialogActions>
                    <Button onClick={this.handleClose} color="primary">
                        Cancel
                    </Button>
                    <Button disabled={name.trim().length === 0} onClick={this.handleSave} color="primary">
                        Save
                    </Button>
                </DialogActions>
            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {
        savedFilter: state.savedDatasetFilter,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleSave: value => {
            dispatch(saveDatasetFilter(value));
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
)(SaveDatasetFilterDialog));


import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';

import DialogTitle from '@mui/material/DialogTitle';
import TextField from '@mui/material/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {saveFeatureSet, setDialog, setMessage} from './actions';

class SaveSetDialog extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {
            open: true,
            name: '',
            category: ''
        };
    }


    onNameChanged = (event) => {
        this.setState({name: event.target.value});
    };

    onCategoryChanged = (event) => {
        this.setState({category: event.target.value});
    };

    handleClose = () => {
        this.props.handleCancel();
    };

    handleSave = () => {
        let name = this.state.name.trim();
        let category = this.state.category;
        if (category != null) {
            category = category.trim();
        }
        this.setState({loading: true});
        this.props.handleSave({name: name, category: category});
    };


    render() {
        let {name, category} = this.state;
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="edit-set-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="edit-set-dialog-title">Create Set</DialogTitle>
                <DialogContent>

                    <TextField
                        size={"small"}
                        required={true}
                        autoComplete="off"
                        value={name}
                        onChange={this.onNameChanged}
                        margin="dense"
                        label="Set Name"
                        fullWidth
                    />
                    <TextField
                        size={"small"}
                        required={true}
                        autoComplete="off"
                        value={category}
                        onChange={this.onCategoryChanged}
                        margin="dense"
                        label="Set Category"
                        fullWidth
                    />
                </DialogContent>
                <DialogActions>
                    <Button onClick={this.handleClose}>
                        Cancel
                    </Button>
                    <Button disabled={name.trim().length === 0 || category.trim().length === 0}
                            onClick={this.handleSave} variant="contained" color="primary">
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
            dispatch(saveFeatureSet(value));
        },
        handleCancel: value => {
            dispatch(setDialog(null));
        },
        handleError: value => {
            dispatch(setMessage(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(SaveSetDialog));


import {CircularProgress} from '@material-ui/core';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {EDIT_DATASET_DIALOG, saveDataset, setDialog, setMessage} from './actions';


function getUniqueArray(text) {
    let tokens = text.split(',');
    let values = new Set();
    for (let i = 0; i < tokens.length; i++) {
        let s = tokens[i].trim().toLowerCase();
        if (s !== '') {
            values.add(s);
        }
    }
    return Array.from(values);
}

class EditDatasetDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            url: '',
            name: this.props.dataset != null ? this.props.dataset.name : '',
            description: this.props.dataset != null ? (this.props.dataset.description != null ? this.props.dataset.description : '') : '',
            readers: '',
            loading: this.props.dataset != null,
        };
    }


    componentDidMount() {

        if (this.props.dataset != null) {
            this.props.serverInfo.api.getDatasetPromise(this.props.dataset.id).then(datasetInfo => {
                let readers = datasetInfo.readers;
                let myIndex = readers.indexOf(this.props.email);
                if (myIndex !== -1) {
                    readers.splice(myIndex, 1);
                }
                this.setState({
                    name: datasetInfo.name,
                    description: this.props.dataset.description != null ? this.props.dataset.description : '',
                    loading: false,
                    url: datasetInfo.url,
                    readers: readers.join(', '),
                });
            }).catch(err => {
                console.log(err);
                this.props.handleError('Unable to retrieve dataset details. Please try again.');
                this.props.handleCancel();
            });
        }
    }

    handleClose = () => {
        this.props.handleCancel();
    };

    handleSave = () => {
        let name = this.state.name.trim();
        let url = this.state.url.trim();
        if (name === '' || url === '') {
            return;
        }
        let description = this.state.description.trim();
        this.setState({loading: true});
        let readers = getUniqueArray(this.state.readers);
        this.props.handleSave({
            dataset: this.props.dataset,
            name: name,
            description: description,
            url: url,
            readers: readers
        });
    };

    onEmailChanged = (event) => {
        this.setState({readers: event.target.value});
    };
    onUrlChanged = (event) => {
        this.setState({url: event.target.value});
    };
    onNameChanged = (event) => {
        this.setState({name: event.target.value});
    };
    onDescriptionChanged = (event) => {
        this.setState({description: event.target.value});
    };

    render() {
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="edit-dataset-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="edit-dataset-dialog-title">{this.props.dataset == null
                    ? 'Import'
                    : 'Edit'} Dataset</DialogTitle>
                <DialogContent>
                    {this.state.loading && <CircularProgress/>}
                    <TextField
                        disabled={this.state.loading}
                        autoComplete="off"
                        required={true}
                        value={this.state.name}
                        onChange={this.onNameChanged}
                        margin="dense"
                        label="Dataset Name"
                        fullWidth
                    />
                    <TextField
                        disabled={this.state.loading}
                        autoComplete="off"
                        required={false}
                        value={this.state.description}
                        onChange={this.onDescriptionChanged}
                        margin="dense"
                        label="Dataset Description"
                        fullWidth
                        multiline={true}
                    />

                    {!this.state.loading &&
                    <TextField
                        required={true}
                        autoComplete="off"
                        value={this.state.url}
                        onChange={this.onUrlChanged}
                        margin="dense"
                        helperText={this.props.serverInfo.serverEmail ? "Please ensure that " + this.props.serverInfo.serverEmail + " has \"Storage Object Viewer\" permissions" : ''}
                        label={"Dataset URL (" + (this.props.serverInfo.serverEmail ? "gs://my_bucket/my_dataset" : "/Users/foo/my_dataset") + ")"}
                        fullWidth
                    />}

                    {!this.state.loading && this.props.serverInfo.clientId !== '' &&
                    <TextField
                        value={this.state.readers}
                        onChange={this.onEmailChanged}
                        margin="dense"
                        label="Share"
                        helperText="Enter comma separated list of emails"
                        fullWidth
                        multiline
                    />}
                </DialogContent>
                <DialogActions>
                    <Button disabled={this.state.loading} onClick={this.handleClose} color="primary">
                        Cancel
                    </Button>
                    <Button disabled={this.state.loading} onClick={this.handleSave} color="primary">
                        Save
                    </Button>
                </DialogActions>


            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dialog === EDIT_DATASET_DIALOG ? state.dataset : null,
        serverInfo: state.serverInfo,
        email: state.email
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleSave: value => {
            dispatch(saveDataset(value));
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
)(EditDatasetDialog));


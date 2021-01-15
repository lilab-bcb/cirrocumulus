import {CircularProgress, InputLabel} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import FormControl from '@material-ui/core/FormControl';
import FormHelperText from '@material-ui/core/FormHelperText';
import Link from '@material-ui/core/Link';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import ReactMarkdown from 'react-markdown';
import {connect} from 'react-redux';
import {EDIT_DATASET_DIALOG, saveDataset, setDialog, setMessage} from './actions';

function TabPanel(props) {
    const {children, value, index, ...other} = props;

    return (
        <div
            role="tabpanel"
            hidden={value !== index}
            id={`editor-tabpanel-${index}`}
            aria-labelledby={`editor-tab-${index}`}
            {...other}
        >
            {value === index && (
                <Box>
                    {children}
                </Box>
            )}
        </div>
    );
}


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
            tabValue: 0,
            name: this.props.dataset != null ? this.props.dataset.name : '',
            title: this.props.title != null ? (this.props.dataset.title != null ? this.props.dataset.title : '') : '',
            description: this.props.dataset != null ? (this.props.dataset.description != null ? this.props.dataset.description : '') : '',
            readers: '',
            loading: this.props.dataset != null,
        };
    }


    componentDidMount() {

        if (this.props.dataset != null) {
            let readers = this.props.dataset.readers;
            if (readers) {
                let myIndex = readers.indexOf(this.props.email);
                if (myIndex !== -1) {
                    readers.splice(myIndex, 1);
                }
                readers = readers.join(', ');
            }
            this.setState({
                name: this.props.dataset.name,
                description: this.props.dataset.description != null ? this.props.dataset.description : '',
                title: this.props.dataset.title != null ? this.props.dataset.title : '',
                loading: false,
                url: this.props.dataset.url,
                readers: readers,
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
        let title = this.state.title.trim();
        this.setState({loading: true});

        let readers = null;
        if (this.state.readers != null) {
            readers = getUniqueArray(this.state.readers);
        }
        this.props.handleSave({
            dataset: this.props.dataset,
            name: name,
            title: title,
            description: description,
            url: url,
            readers: readers
        });
    };

    onTabChanged = (event, value) => {
        this.setState({tabValue: value});
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
    onTitleChanged = (event) => {
        this.setState({title: event.target.value});
    };

    render() {
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="edit-dataset-dialog-title"
                fullWidth={true}
                maxWidth={'lg'}
            >
                <DialogTitle id="edit-dataset-dialog-title">{this.props.dataset == null
                    ? 'New'
                    : 'Edit'} Dataset</DialogTitle>
                <DialogContent>
                    {this.state.loading && <CircularProgress/>}
                    <TextField
                        disabled={this.state.loading}
                        autoComplete="off"
                        required={true}
                        value={this.state.name}
                        helperText={"asfd"}
                        onChange={this.onNameChanged}
                        margin="dense"
                        label="Name"
                        fullWidth
                    />

                    {!this.state.loading &&
                    <TextField
                        required={true}
                        autoComplete="off"
                        value={this.state.url}
                        onChange={this.onUrlChanged}
                        margin="dense"
                        helperText={this.props.serverInfo.email ? "Please ensure that " + this.props.serverInfo.email + " has \"Storage Object Viewer\" permissions" : ''}
                        label={"URL (" + (this.props.serverInfo.email ? "gs://my_bucket/my_dataset" : "/Users/foo/my_dataset") + ")"}
                        fullWidth
                    />}

                    <TextField
                        disabled={this.state.loading}
                        autoComplete="off"
                        required={false}
                        value={this.state.title}
                        onChange={this.onTitleChanged}
                        margin="dense"
                        label="Title"
                        fullWidth
                        inputProps={{maxLength: 255}}
                    />
                    <div>
                        <FormControl>
                            <InputLabel>Summary</InputLabel>
                            <div style={{marginTop: 36}}></div>
                            <FormHelperText><Link
                                href={"https://www.markdownguide.org/cheat-sheet/"}
                                target="_blank">Markdown Cheat Sheet</Link></FormHelperText>
                        </FormControl>

                    </div>

                    <Tabs value={this.state.tabValue} onChange={this.onTabChanged}
                          aria-label="description-editor">
                        <Tab label="Write"/>
                        <Tab label="Preview"/>
                    </Tabs>

                    <TabPanel value={this.state.tabValue} index={0}>
                        <TextField
                            disabled={this.state.loading}
                            autoComplete="off"
                            required={false}
                            value={this.state.description}
                            onChange={this.onDescriptionChanged}
                            margin="dense"
                            fullWidth
                            variant="outlined"
                            rows={8}
                            rowsMax={8}
                            multiline={true}
                            inputProps={{maxLength: 1000}}
                        />
                    </TabPanel>
                    <TabPanel value={this.state.tabValue} index={1}>
                        {this.state.description !== '' && <Box border={1}>
                            <ReactMarkdown linkTarget="_blank" children={this.state.description}/>
                        </Box>}
                    </TabPanel>


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
        email: state.email,
        serverInfo: state.serverInfo,
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


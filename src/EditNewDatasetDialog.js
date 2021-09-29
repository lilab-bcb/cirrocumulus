import {CircularProgress, InputLabel, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';

import DialogTitle from '@mui/material/DialogTitle';
import Divider from '@mui/material/Divider';
import FormControl from '@mui/material/FormControl';
import FormHelperText from '@mui/material/FormHelperText';
import Link from '@mui/material/Link';
import MenuItem from '@mui/material/MenuItem';
import Select from '@mui/material/Select';
import withStyles from '@mui/styles/withStyles';
import Tab from '@mui/material/Tab';
import Tabs from '@mui/material/Tabs';
import TextField from '@mui/material/TextField';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';
import LinkIcon from '@mui/icons-material/Link';
import ReactMarkdown from 'markdown-to-jsx';
import React from 'react';
import {connect} from 'react-redux';
import {EDIT_DATASET_DIALOG, saveDataset, setDialog, setMessage} from './actions';
import {REACT_MD_OVERRIDES} from './util';

const styles = theme => ({

    formControl: {
        minWidth: 200,
        marginTop: theme.spacing(1)
    },
    select: {
        minWidth: 200
    }
});

const favoriteSpecies = ["Homo sapiens", "Mus musculus"];
const otherSpecies = ["Gallus gallus", "Macaca fascicularis", "Macaca mulatta", "Rattus norvegicus"];


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

class EditNewDatasetDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.init = false;
        this.fileInputRef = React.createRef();
        this.fileDropRef = React.createRef();
        this.dragIndicator = React.createRef();
        this.state = {
            url: '',
            writePreviewTabValue: "Write",
            uploadTabValue: "URL",
            species: this.props.dataset != null ? this.props.dataset.species : '',
            name: this.props.dataset != null ? this.props.dataset.name : '',
            title: this.props.title != null ? (this.props.dataset.title != null ? this.props.dataset.title : '') : '',
            description: this.props.dataset != null ? (this.props.dataset.description != null ? this.props.dataset.description : '') : '',
            readers: '',
            loading: this.props.dataset != null
        };
    }

    showDragIndicator(event, show) {
        event.stopPropagation();
        event.preventDefault();
        this.dragIndicator.current.style.background = show ? '#1976d2' : '';
    }


    initDragDrop = () => {
        if (!this.init) {
            const fileDropRef = this.fileDropRef.current;
            this.init = true;

            fileDropRef.addEventListener('dragover', (e) => {
                this.showDragIndicator(e, true);
            });
            fileDropRef.addEventListener('dragenter', (e) => {
                this.showDragIndicator(e, true);
            });
            fileDropRef.addEventListener('dragend', (e) => {
                this.showDragIndicator(e, false);
            });
            fileDropRef.addEventListener('dragleave', (e) => {
                this.showDragIndicator(e, false);
            });
            fileDropRef.addEventListener('drop', (e) => {
                this.showDragIndicator(e, false);
                const filesArray = e.dataTransfer.files;
                if (filesArray.length === 1) {
                    this.setState({file: filesArray[0]});
                }
            });
        }
    };


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
                file: null,
                species: this.props.dataset.species != null ? this.props.dataset.species : '',
                name: this.props.dataset.name,
                description: this.props.dataset.description != null ? this.props.dataset.description : '',
                title: this.props.dataset.title != null ? this.props.dataset.title : '',
                loading: false,
                url: this.props.dataset.url,
                readers: readers
            });
        }
    }

    handleClose = () => {
        this.props.handleCancel();
    };


    handleSave = () => {
        let name = this.state.name.trim();
        if (name === '') {
            return;
        }
        let url = null;
        let file = null;
        const isNew = this.props.dataset == null;
        if (isNew) {
            if (this.state.uploadTabValue === 1 || !this.props.serverInfo.upload) {
                url = this.state.url.trim();

            } else {
                file = this.state.file;
            }
        }

        let description = this.state.description.trim();
        let title = this.state.title.trim();
        let species = this.state.species;
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
            file: file,
            species: species,
            readers: readers
        });
    };

    onWritePreviewTabChanged = (event, value) => {
        this.setState({writePreviewTabValue: value});
    };

    onUploadTabChanged = (event, value) => {
        this.setState({uploadTabValue: value});
    };

    onEmailChanged = (event) => {
        this.setState({readers: event.target.value});
    };
    onUrlChanged = (event) => {
        this.setState({url: event.target.value});
    };
    onFilesChanged = () => {
        const files = this.fileInputRef.current.files;
        this.setState({file: files[0]});
    };
    onSpeciesChange = (event) => {
        this.setState({species: event.target.value});
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
        const canUpload = this.props.serverInfo.upload;
        const isNew = this.props.dataset == null;
        let saveEnabled = !this.state.loading && this.state.name.trim() !== '';
        const isAuthEnabled = this.props.serverInfo.clientId !== '';

        if (isNew) {
            if (this.state.uploadTabValue === 1 || !this.props.serverInfo.upload) {
                const url = this.state.url.trim();
                saveEnabled = saveEnabled && url !== '';
            } else {
                saveEnabled = saveEnabled && this.state.file != null;
            }
        }
        return (
            <Dialog
                open={true}
                ref={this.fileDropRef}
                onClose={this.handleClose}
                aria-labelledby="edit-dataset-dialog-title"
                fullWidth={true}
                maxWidth={'lg'}
                TransitionProps={{
                    onEntered: e => this.initDragDrop()
                }}>
                <DialogTitle id="edit-dataset-dialog-title">{isNew
                    ? 'New'
                    : 'Edit'} Dataset</DialogTitle>
                <DialogContent>
                    <Box>
                        {this.state.loading && <CircularProgress/>}
                        <TextField
                            size={"small"}
                            disabled={this.state.loading}
                            autoComplete="off"
                            required={true}
                            value={this.state.name}
                            onChange={this.onNameChanged}
                            margin="dense"
                            label="Name"
                            fullWidth
                        />

                        <div style={{display: isNew && canUpload ? '' : 'none'}}>
                            <FormControl className={this.props.classes.formControl}>
                                <InputLabel shrink required>Source</InputLabel>
                            </FormControl>
                            <Tabs value={this.state.uploadTabValue} onChange={this.onUploadTabChanged}
                                  aria-label="upload">
                                <Tab value="My Computer" label="My Computer" icon={<CloudUploadIcon/>}/>
                                <Tab value="URL" label="URL" icon={<LinkIcon/>}/>
                            </Tabs>
                            <div
                                role="tabpanel"
                                hidden={this.state.uploadTabValue !== "My Computer"}
                            >
                                <div ref={this.dragIndicator}>
                                    <Button size="small" variant="outlined" disabled={this.state.loading}
                                            onClick={e => this.fileInputRef.current.click()}>Select File</Button>
                                    <Typography style={{display: 'inline-block', paddingLeft: '1em'}} component={"h3"}>or
                                        Drag
                                        And Drop</Typography>
                                    <input hidden ref={this.fileInputRef} type="file" onChange={this.onFilesChanged}/>
                                    <Typography style={{display: 'block'}} color="textPrimary"
                                                variant={"caption"}>{this.state.file ? this.state.file.name : ''}</Typography>
                                </div>
                                <Divider style={{marginTop: '1em', marginBottom: '1em'}}/>
                            </div>


                            <div role="tabpanel" hidden={this.state.uploadTabValue !== "URL"}>
                                <TextField
                                    size={"small"}
                                    required={true}
                                    disabled={this.state.loading}
                                    autoComplete="off"
                                    value={this.state.url}
                                    onChange={this.onUrlChanged}
                                    margin="dense"
                                    helperText={this.props.serverInfo.email ? "Please ensure that " + this.props.serverInfo.email + " has read permission to this URL" : ''}
                                    label={"URL"}
                                    fullWidth
                                />
                            </div>
                        </div>

                        <TextField
                            size={"small"}
                            style={{display: isNew && canUpload ? 'none' : ''}}
                            required={true}
                            disabled={this.state.loading || !isNew}
                            autoComplete="off"
                            value={this.state.url}
                            onChange={this.onUrlChanged}
                            margin="dense"
                            helperText={this.props.serverInfo.email ? "Please ensure that " + this.props.serverInfo.email + " has read permission to this URL" : ''}
                            label={"URL"}
                            fullWidth
                        />

                        <FormControl className={this.props.classes.formControl}>
                            <InputLabel id="species-label">Species</InputLabel>
                            <Select
                                label={"Species"}
                                size={"small"}
                                labelId="species-label"
                                value={this.state.species}
                                onChange={this.onSpeciesChange}
                            >
                                {favoriteSpecies.map(species => <MenuItem key={species}
                                                                          value={species}>{species}</MenuItem>)}
                                <MenuItem divider={true}/>
                                {otherSpecies.map(species => <MenuItem key={species}
                                                                       value={species}>{species}</MenuItem>)}
                            </Select>
                        </FormControl>


                        <TextField
                            size={"small"}
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
                        <FormControl className={this.props.classes.formControl}>
                            <InputLabel shrink>Summary</InputLabel>
                            <FormHelperText style={{marginTop: 14}}><Link
                                href={"https://www.markdownguide.org/cheat-sheet/"}
                                target="_blank">Markdown Cheat Sheet</Link></FormHelperText>
                        </FormControl>


                        <Tabs value={this.state.writePreviewTabValue} onChange={this.onWritePreviewTabChanged}
                              aria-label="description-editor">
                            <Tab value="Write" label="Write"/>
                            <Tab value="Preview" label="Preview"/>
                        </Tabs>

                        <div role="tabpanel" hidden={this.state.writePreviewTabValue !== "Write"}>
                            <TextField
                                size={"small"}
                                disabled={this.state.loading}
                                autoComplete="off"
                                required={false}
                                value={this.state.description}
                                onChange={this.onDescriptionChanged}
                                margin="dense"
                                fullWidth
                                variant="outlined"
                                rows={8}
                                maxRows={8}
                                multiline={true}
                                inputProps={{maxLength: 1000}}
                            />
                        </div>
                        <div role="tabpanel" hidden={this.state.writePreviewTabValue !== "Preview"}>
                            {this.state.description !== '' && <Box border={1}>
                                <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                               children={this.state.description}/>
                            </Box>}
                        </div>

                        <TextField
                            size={"small"}
                            value={this.state.readers}
                            onChange={this.onEmailChanged}
                            margin="dense"
                            label="Readers"
                            helperText={isAuthEnabled ? "Enter comma separated list of emails" : "Enable authentication to permit sharing"}
                            disabled={this.state.loading || !isAuthEnabled}
                            fullWidth
                            multiline
                        />
                    </Box>
                </DialogContent>
                <DialogActions>
                    <Button disabled={this.state.loading} onClick={this.handleClose}>
                        Cancel
                    </Button>
                    <Button disabled={!saveEnabled} onClick={this.handleSave}
                            variant="contained" color="primary">
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
        serverInfo: state.serverInfo
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
        }

    };
};


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(EditNewDatasetDialog));

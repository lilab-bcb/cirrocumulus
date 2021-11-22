import {InputLabel, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import DialogActions from '@mui/material/DialogActions';
import DialogContent from '@mui/material/DialogContent';
import LoadingButton from '@mui/lab/LoadingButton';
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
import React, {useCallback, useEffect, useRef, useState} from 'react';
import {connect} from 'react-redux';
import {EDIT_DATASET_DIALOG, saveDataset, setDialog} from './actions';
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

export const DATASET_FIELDS = [{fieldName: 'title', label: 'Title'}, {
    fieldName: 'species',
    label: 'Species'
}, {fieldName: 'experimentType', label: 'Experiment Type'}, {fieldName: 'description', label: 'Summary'},
    {fieldName: 'overallDesign', label: 'Overall Design'}, {fieldName: 'citations', label: 'Citation(s)'}];

function EditNewDatasetDialog(props) {
    const {classes, dataset, email, serverInfo, onSave, onCancel} = props;

    const fileInputRef = useRef();
    const fileDropRef = useRef();
    const dragIndicator = useRef();

    const [file, setFile] = useState(null);
    const [writePreviewTabValue, setWritePreviewTabValue] = useState('Write');
    const [uploadTabValue, setUploadTabValue] = useState('URL');
    const [url, setUrl] = useState('');
    const [loading, setLoading] = useState(false);
    const [readers, setReaders] = useState('');

    const [name, setName] = useState('');
    const [title, setTitle] = useState('');
    const [species, setSpecies] = useState('');
    const [experimentType, setExperimentType] = useState('');
    const [description, setDescription] = useState(''); // summary
    const [overallDesign, setOverallDesign] = useState('');
    const [citations, setCitations] = useState('');

    const favoriteSpecies = serverInfo.species.favorite;
    const otherSpecies = serverInfo.species.other;


    const canUpload = serverInfo.upload;
    const isNew = dataset == null;
    let saveEnabled = !loading && name.trim() !== '';
    const isAuthEnabled = serverInfo.clientId !== '';

    if (isNew) {
        if (uploadTabValue === 'URL' || !serverInfo.upload) {
            saveEnabled = saveEnabled && url.trim() !== '';
        } else {
            saveEnabled = saveEnabled && file != null;
        }
    }

    const setFileDropRef = useCallback(node => {

        function dragStart(event) {
            showDragIndicator(event, true);
        }

        function dragEnd(event) {
            showDragIndicator(event, false);
        }

        function drop(event) {
            showDragIndicator(event, false);
            const filesArray = event.dataTransfer.files;
            if (filesArray.length === 1) {
                setFile(filesArray[0]);
            }
        }

        function showDragIndicator(event, show) {
            event.stopPropagation();
            event.preventDefault();
            dragIndicator.current.style.background = show ? '#1976d2' : '';
        }


        if (canUpload && isNew) {
            if (fileDropRef.current) {
                fileDropRef.current.removeEventListener('dragover', dragStart);
                fileDropRef.current.removeEventListener('dragenter', dragStart);
                fileDropRef.current.removeEventListener('dragend', dragEnd);
                fileDropRef.current.removeEventListener('dragleave', dragEnd);
                fileDropRef.current.removeEventListener('drop', drop);
            }
            fileDropRef.current = node;
            if (node) {
                fileDropRef.current.addEventListener('dragover', dragStart);
                fileDropRef.current.addEventListener('dragenter', dragStart);
                fileDropRef.current.addEventListener('dragend', dragEnd);
                fileDropRef.current.addEventListener('dragleave', dragEnd);
                fileDropRef.current.addEventListener('drop', drop);
            }
        }

    }, [canUpload, isNew]);


    useEffect(() => {
        if (dataset != null) {
            let readers = dataset.readers;
            if (readers) {
                let myIndex = readers.indexOf(email);
                if (myIndex !== -1) {
                    readers.splice(myIndex, 1);
                }
                readers = readers.join(', ');
            }
            setReaders(readers);
        } else {
            setReaders('');
        }
        setLoading(false);
        setUrl(dataset != null ? dataset.url : '');


        setName(dataset != null && dataset.name != null ? dataset.name : '');
        setTitle(dataset != null && dataset.title != null ? dataset.title : '');
        setSpecies(dataset != null && dataset.species != null ? dataset.species : '');
        setExperimentType(dataset != null && dataset.experimentType != null ? dataset.experimentType : '');
        setDescription(dataset != null && dataset.description != null ? dataset.description : '');
        setOverallDesign(dataset != null && dataset.overallDesign != null ? dataset.overallDesign : '');
        setCitations(dataset != null && dataset.citations != null ? dataset.citations : '');

    }, [dataset, email]);


    function handleSave() {
        if (name.trim() === '') {
            return;
        }
        let useUrl = false;
        const isNew = dataset == null;
        if (isNew) {
            if (uploadTabValue === 'URL' || !serverInfo.upload) {
                useUrl = true;
            }
        }

        setLoading(true);
        onSave({
            dataset: dataset,
            url: isNew && useUrl ? url : null,
            file: isNew && !useUrl ? file : null,
            readers: readers.trim() !== '' ? getUniqueArray(readers) : null,
            name: name.trim(),
            title: title.trim(),
            species: species,
            experimentType: experimentType.trim(),
            description: description.trim(),
            overallDesign: overallDesign.trim(),
            citations: citations.trim()
        });
    }


    return (
        <Dialog
            open={true}
            ref={setFileDropRef}
            onClose={onCancel}
            aria-labelledby="edit-dataset-dialog-title"
            fullWidth={true}
            maxWidth={'lg'}>
            <DialogTitle id="edit-dataset-dialog-title">{isNew
                ? 'New'
                : 'Edit'} Dataset</DialogTitle>
            <DialogContent>
                <Box>
                    <TextField
                        size={"small"}
                        disabled={loading}
                        autoComplete="off"
                        required={true}
                        value={name}
                        onChange={event => setName(event.target.value)}
                        margin="dense"
                        label="Name"
                        fullWidth
                    />

                    <div style={{display: isNew && canUpload ? '' : 'none'}}>
                        <FormControl className={classes.formControl}>
                            <InputLabel shrink required>Source</InputLabel>
                        </FormControl>
                        <Tabs disabled={loading} value={uploadTabValue}
                              onChange={(event, value) => setUploadTabValue(value)}
                              aria-label="upload">
                            <Tab value="My Computer" label="My Computer" icon={<CloudUploadIcon/>}/>
                            <Tab value="URL" label="URL" icon={<LinkIcon/>}/>
                        </Tabs>
                        <div
                            role="tabpanel"
                            hidden={uploadTabValue !== "My Computer"}
                        >
                            <div ref={dragIndicator}>
                                <Button size="small" variant="outlined" disabled={loading}
                                        onClick={e => fileInputRef.current.click()}>Select File</Button>
                                <Typography style={{display: 'inline-block', paddingLeft: '1em'}} component={"h3"}>or
                                    Drag
                                    And Drop</Typography>
                                <input hidden ref={fileInputRef} type="file"
                                       onChange={event => setFile(fileInputRef.current.files[0])}/>
                                <Typography style={{display: 'block'}} color="textPrimary"
                                            variant={"caption"}>{file ? file.name : ''}</Typography>
                            </div>
                            <Divider style={{marginTop: '1em', marginBottom: '1em'}}/>
                        </div>


                        <div role="tabpanel" hidden={uploadTabValue !== "URL"}>
                            <TextField
                                size={"small"}
                                required={true}
                                disabled={loading}
                                autoComplete="off"
                                value={url}
                                onChange={event => setUrl(event.target.value)}
                                margin="dense"
                                helperText={serverInfo.email ? "Please ensure that " + serverInfo.email + " has read permission to this URL" : ''}
                                label={"URL"}
                                fullWidth
                            />
                        </div>
                    </div>
                    <TextField
                        size={"small"}
                        style={{display: isNew && canUpload ? 'none' : ''}}
                        required={true}
                        disabled={loading || !isNew}
                        autoComplete="off"
                        value={url}
                        onChange={event => setUrl(event.target.value)}
                        margin="dense"
                        helperText={serverInfo.email ? "Please ensure that " + serverInfo.email + " has read permission to this URL" : ''}
                        label={"URL"}
                        fullWidth
                    />

                    <TextField
                        size={"small"}
                        disabled={loading}
                        autoComplete="off"
                        required={false}
                        value={title}
                        onChange={event => setTitle(event.target.value)}
                        margin="dense"
                        label="Title"
                        fullWidth
                        inputProps={{maxLength: 255}}
                    />
                    <FormControl className={classes.formControl}>
                        <InputLabel id="species-label">Species</InputLabel>
                        <Select
                            disabled={loading}
                            label={"Species"}
                            size={"small"}
                            labelId="species-label"
                            value={species}
                            onChange={event => setSpecies(event.target.value)}
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
                        disabled={loading}
                        autoComplete="off"
                        value={experimentType}
                        onChange={(event) => setExperimentType(event.target.value)}
                        margin="dense"
                        label={"Experiment Type"}
                        fullWidth
                    />
                    <FormControl className={classes.formControl}>
                        <InputLabel shrink>Summary</InputLabel>
                        <FormHelperText style={{marginTop: 14}}><Link
                            href={"https://www.markdownguide.org/cheat-sheet/"}
                            target="_blank">Markdown Cheat Sheet</Link></FormHelperText>
                    </FormControl>
                    <Tabs disabled={loading} value={writePreviewTabValue}
                          onChange={(event, value) => setWritePreviewTabValue(value)}
                          aria-label="description-editor">
                        <Tab value="Write" label="Write"/>
                        <Tab value="Preview" label="Preview"/>
                    </Tabs>
                    <div role="tabpanel" hidden={writePreviewTabValue !== "Write"}>
                        <TextField
                            size={"small"}
                            disabled={loading}
                            autoComplete="off"
                            required={false}
                            value={description}
                            onChange={event => setDescription(event.target.value)}
                            margin="dense"
                            fullWidth
                            variant="outlined"
                            rows={8}
                            maxRows={8}
                            multiline={true}
                            inputProps={{maxLength: 1000}}
                        />
                    </div>
                    <div role="tabpanel" hidden={writePreviewTabValue !== "Preview"}>
                        {description !== '' && <Box border={1}>
                            <ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}}
                                           children={description}/>
                        </Box>}
                    </div>
                    <TextField
                        size={"small"}
                        disabled={loading}
                        autoComplete="off"
                        required={false}
                        value={overallDesign}
                        onChange={(event) => setOverallDesign(event.target.value)}
                        margin="dense"
                        label="Overall Design"
                        fullWidth
                    />
                    <TextField
                        size={"small"}
                        disabled={loading}
                        autoComplete="off"
                        value={citations}
                        onChange={(event) => setCitations(event.target.value)}
                        margin="dense"
                        label={"Citation(s)"}
                        fullWidth
                    />


                    <TextField
                        size={"small"}
                        value={readers}
                        onChange={event => setReaders(event.target.value)}
                        margin="dense"
                        label="Readers"
                        helperText={isAuthEnabled ? "Enter comma separated list of emails" : "Enable authentication to permit sharing"}
                        disabled={loading || !isAuthEnabled}
                        fullWidth
                        multiline
                    />
                </Box>
            </DialogContent>
            <DialogActions>
                <Button disabled={loading} onClick={onCancel}>
                    Cancel
                </Button>
                <LoadingButton disabled={!saveEnabled} onClick={handleSave}
                               variant="contained" color="primary" loading={loading}>
                    Save
                </LoadingButton>
            </DialogActions>
        </Dialog>
    );

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
        onSave: value => {
            dispatch(saveDataset(value));
        },
        onCancel: value => {
            dispatch(setDialog(null));
        }
    };
};


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps
)(EditNewDatasetDialog));

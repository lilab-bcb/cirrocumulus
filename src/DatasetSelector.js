import {Typography} from '@material-ui/core';
import Button from '@material-ui/core/Button';
import CircularProgress from '@material-ui/core/CircularProgress';
import IconButton from '@material-ui/core/IconButton';
import InputAdornment from '@material-ui/core/InputAdornment';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import Popover from '@material-ui/core/Popover';
import TextField from '@material-ui/core/TextField';
import ClearIcon from '@material-ui/icons/Clear';
import FolderOpenIcon from '@material-ui/icons/FolderOpen';
import InfoIcon from '@material-ui/icons/Info';
import {groupBy} from 'lodash';
import React from 'react';
import ReactMarkdown from 'react-markdown';
import {connect} from 'react-redux';
import {setDialog} from './actions';

export class DatasetSelector extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {anchorEl: null, searchText: '', datasetDetailsEl: null, selectedDataset: null};
    }

    handleClick = (event) => {
        this.setState({anchorEl: event.currentTarget, searchText: ''});
    };

    handleClose = () => {
        this.setState({anchorEl: null, searchText: ''});
    };

    handleCloseDatasetDetails = (event) => {
        this.setState({datasetDetailsEl: null, selectedDataset: null});
    };

    handleListItemDetailsClick = (event, id) => {
        this.setState({datasetDetailsEl: event.currentTarget});
        this.props.serverInfo.api.getDatasetPromise(id).then(result => {
            this.setState({selectedDataset: result});
        });

    };

    handleListItemClick = (id) => {
        const selectedId = this.props.dataset != null ? this.props.dataset.id : null;
        if (id !== selectedId) {
            this.props.onChange(id);
        }
        this.setState({anchorEl: null, searchText: ''});
    };

    handleClearSearchText = () => {
        this.setState({searchText: ''});
    };

    onSearchChange = (event) => {
        this.setState({searchText: event.target.value});
    };

    render() {
        const {dataset, datasetChoices} = this.props;
        const selectedId = dataset != null ? dataset.id : null;

        if (datasetChoices.length <= 1 && selectedId != null) {
            return null;
        }
        const {anchorEl, searchText, selectedDataset} = this.state;

        const open = Boolean(this.state.anchorEl);
        let filteredChoices = datasetChoices;
        const searchTextLower = this.state.searchText.toLowerCase().trim();
        if (searchTextLower != '') {
            filteredChoices = filteredChoices.filter(choice => choice.name.toLowerCase().indexOf(searchTextLower) !== -1 ||
                (choice.description != null && choice.description.toLowerCase().indexOf(searchTextLower) !== -1));
        }
        const datasetDetailsOpen = Boolean(this.state.datasetDetailsEl);
        const hasMoreInfo = selectedDataset && (selectedDataset.title || selectedDataset.description);
        const species2Items = groupBy(filteredChoices, 'species');
        const speciesArray = Object.keys(species2Items);
        speciesArray.sort();
        return (
            <React.Fragment>
                <Popover
                    id={"dataset-details-selector"}
                    open={datasetDetailsOpen}
                    anchorEl={this.state.datasetDetailsEl}
                    onClose={this.handleCloseDatasetDetails}
                    anchorOrigin={{
                        vertical: 'bottom',
                        horizontal: 'center',
                    }}
                    transformOrigin={{
                        vertical: 'top',
                        horizontal: 'center',
                    }}
                >
                    <div style={{width: 500}}>
                        {selectedDataset == null &&
                        <div><CircularProgress size={20}/> Loading...</div>}
                        {selectedDataset && <Typography component={"h6"}>
                            {selectedDataset.name}
                        </Typography>}
                        {!hasMoreInfo && <div>No description available</div>}

                        {selectedDataset && selectedDataset.title && <div>{selectedDataset.title}</div>}
                        {selectedDataset && selectedDataset.description &&
                        <ReactMarkdown linkTarget="_blank" children={selectedDataset.description}/>}
                    </div>
                </Popover>
                <Button variant="contained" onClick={this.handleClick}
                        color={selectedId == null ? "primary" : "default"}
                        startIcon={<FolderOpenIcon/>}>Open</Button>
                <Popover open={open}
                         anchorEl={anchorEl}
                         onClose={this.handleClose}
                         anchorOrigin={{
                             vertical: 'bottom',
                             horizontal: 'center',
                         }}
                         transformOrigin={{
                             vertical: 'top',
                             horizontal: 'center',
                         }}
                >
                    <TextField style={{padding: 6}} type="text" placeholder={"Search"} value={searchText}
                               onChange={this.onSearchChange}
                               fullWidth={true}
                               InputProps={searchText.trim() !== '' ? {
                                   endAdornment:
                                       <InputAdornment position="end">
                                           <IconButton
                                               aria-label="clear"
                                               onClick={this.handleClearSearchText}
                                           >
                                               <ClearIcon/>
                                           </IconButton>
                                       </InputAdornment>
                               } : null}
                    />
                    {speciesArray.map(species => {
                        const speciesText = species === 'null' ? 'Other' : species;
                        const choices = species2Items[species];
                        return <React.Fragment key={species}>
                            <Typography component={"h2"}>{speciesText}</Typography>
                            <List style={{width: 500}} dense disablePadding component="nav" aria-label="datasets">
                                {choices.map(choice => {
                                    let text = choice.name;
                                    if (choice.title) {
                                        text += ' - ' + choice.title;
                                    }
                                    return <ListItem alignItems="flex-start" selected={choice.id === selectedId}
                                                     key={choice.id}
                                                     button onClick={(e) => this.handleListItemClick(choice.id)}>
                                        <ListItemText
                                            primary={text}
                                            style={{
                                                textOverflow: 'ellipsis',
                                                overflow: 'hidden',
                                                whiteSpace: 'nowrap'
                                            }}/>
                                        <ListItemSecondaryAction>
                                            <IconButton onClick={(e) => this.handleListItemDetailsClick(e, choice.id)}
                                                        edge="end"
                                                        aria-label="summary">
                                                <InfoIcon/>
                                            </IconButton>
                                        </ListItemSecondaryAction>
                                    </ListItem>;
                                })}
                            </List></React.Fragment>;
                    })}
                </Popover>
            </React.Fragment>
        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
        serverInfo: state.serverInfo
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleClose: value => {
            dispatch(setDialog(null));
        }

    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DatasetSelector));


import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import InputAdornment from '@material-ui/core/InputAdornment';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import Popover from '@material-ui/core/Popover';
import TextField from '@material-ui/core/TextField';
import ClearIcon from '@material-ui/icons/Clear';
import FolderOpenIcon from '@material-ui/icons/FolderOpen';
import React from 'react';
import {connect} from 'react-redux';
import {setDialog} from './actions';

export class DatasetSelector extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {anchorEl: null, searchText: ''};
    }

    handleClick = (event) => {
        this.setState({anchorEl: event.currentTarget, searchText: ''});
    };

    handleClose = () => {
        this.setState({anchorEl: null, searchText: ''});
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
        const {anchorEl, searchText} = this.state;

        const selectedId = dataset != null ? dataset.id : null;
        const open = Boolean(this.state.anchorEl);
        let filteredChoices = datasetChoices;
        const searchTextLower = this.state.searchText.toLowerCase().trim();
        if (searchTextLower != '') {
            filteredChoices = filteredChoices.filter(choice => choice.name.toLowerCase().indexOf(searchTextLower) !== -1 ||
                (choice.description != null && choice.description.toLowerCase().indexOf(searchTextLower) !== -1));
        }

        return (
            <React.Fragment>
                <Button variant="contained" onClick={this.handleClick}
                        color={selectedId == null ? "primary" : "inherit"}
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
                    <List style={{width: 500}} dense disablePadding component="nav" aria-label="datasets">
                        {filteredChoices.map(choice => {
                            return <ListItem alignItems="flex-start" selected={choice.id === selectedId} key={choice.id}
                                             button onClick={(e) => this.handleListItemClick(choice.id)}>
                                <ListItemText
                                    primary={choice.name} secondary={choice.description}/>
                            </ListItem>;
                        })}
                    </List>
                </Popover>
            </React.Fragment>
        );
    }
}

const mapStateToProps = state => {
    return {
        dataset: state.dataset,
        datasetChoices: state.datasetChoices,
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


import Button from '@material-ui/core/Button';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import Popover from '@material-ui/core/Popover';
import FolderOpenIcon from '@material-ui/icons/FolderOpen';
import React from 'react';
import {connect} from 'react-redux';
import {setDialog} from './actions';

export class DatasetSelector extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {anchorEl: null};
    }

    handleClick = (event) => {
        this.setState({anchorEl: event.currentTarget});
    };

    handleClose = () => {
        this.setState({anchorEl: null});
    };

    handleListItemClick = (id) => {
        const selectedId = this.props.dataset != null ? this.props.dataset.id : null;
        if (id !== selectedId) {
            this.props.onChange(id);
        }
        this.setState({anchorEl: null});
    };

    render() {
        const {datasetChoices, dataset} = this.props;
        const {anchorEl} = this.state;
        const selectedId = dataset != null ? dataset.id : null;
        const open = Boolean(this.state.anchorEl);
        if (datasetChoices.length <= 1 && dataset != null) {
            return null;
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
                    <List style={{maxWidth: 500}} dense disablePadding component="nav" aria-label="datasets">
                        {datasetChoices.map(choice => {
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


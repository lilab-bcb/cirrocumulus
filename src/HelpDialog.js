import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import Divider from '@material-ui/core/Divider';
import Typography from '@material-ui/core/Typography';
import React from 'react';
import {connect} from 'react-redux';
import {setDialog} from './actions';

class HelpDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            open: true,
        };
    }

    handleClose = () => {
        this.props.handleClose();
    };


    render() {
        return (
            <Dialog
                open={this.state.open}
                onClose={this.handleClose}
                aria-labelledby="help-dataset-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="help-dataset-dialog-title">Help</DialogTitle>
                <DialogContent>
                    <Typography variant="h6">Primary View Controls</Typography>
                    <Typography variant="subtitle1">3D</Typography>
                    Rotate: Mouse left click.<br/>
                    Pan: Mouse right click.<br/>
                    Zoom: Mouse wheel.<br/>
                    Holding ctrl reverses the mouse clicks.<br/>
                    <Typography variant="subtitle1">2D</Typography>
                    Pan: Mouse left click.<br/>
                    Zoom: Mouse wheel.<br/>
                    <Divider/>
                    <Typography variant="h6">Gallery</Typography>
                    Drag and drop to reorder the plots. Click the id link to make it the primary view.
                    <Divider/>
                    <Typography variant="h6">Dot Plot</Typography>
                    Enter one or more categorical observations and one or more features.
                </DialogContent>

            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {};
};
const mapDispatchToProps = dispatch => {
    return {


        handleClose: value => {
            dispatch(setDialog(null));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(HelpDialog));


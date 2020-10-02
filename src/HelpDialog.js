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
                    <Typography variant="h6">Primary Embedding Controls</Typography>
                    <Typography variant="subtitle1">3-d Plot</Typography>
                    Rotate: Mouse left click.<br/>
                    Pan: Mouse right click.<br/>
                    Zoom: Mouse wheel.<br/>

                    <Typography variant="subtitle1">2-d Plot</Typography>
                    Pan: Mouse left click.<br/>
                    Zoom: Mouse wheel.<br/>
                    <Divider/>
                    <Typography variant="h6">Gallery</Typography>
                    Drag chart to reorder. Click plot to set primary view.
                    <Divider/>
                    <Typography variant="h6">Dot Plot</Typography>
                    Enter one or more categorical observations and one or more features.
                    <Divider/>
                    <Typography variant="h6">Additional documentation is available <a target="_blank"
                                                                                      rel="noopener noreferrer"
                                                                                      href="http://cirrocumulus.readthedocs.io/">here</a>.</Typography>
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


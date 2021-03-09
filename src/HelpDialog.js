import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import React from 'react';
import {connect} from 'react-redux';
import {setDialog} from './actions';
import LandingPage from './LandingPage';

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
                <DialogContent>
                    <LandingPage/>
                    {/*<Divider/>*/}

                    {/*<small>Last Updated {CIRRO_BUILD_DATE}</small>*/}
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


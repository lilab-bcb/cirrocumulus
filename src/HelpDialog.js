import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import Link from '@material-ui/core/Link';
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
                <DialogContent>
                    <Typography variant="h4">Links</Typography>

                    <div style={{marginLeft: 6, marginBottom: 12}}>
                        <Typography variant="h6"><Link target="_blank" rel="noopener noreferrer"
                                                       href="http://cirrocumulus.readthedocs.io/">Documentation</Link></Typography>
                        <Typography variant="h6">Contact: cirrocumulus@broadinstitute.org</Typography>
                        <Typography variant="h6"><Link target="_blank" rel="noopener noreferrer"
                                                       href="https://github.com/klarman-cell-observatory/cirrocumulus">Source
                            Code</Link></Typography>
                    </div>

                    <Typography variant="h4">Primary Embedding Controls</Typography>
                    <div style={{marginLeft: 6, marginBottom: 12}}>
                        <Typography variant="h5">3-d Embedding</Typography>
                        Rotate: Mouse left click<br/>
                        Pan: Mouse right click<br/>
                        Zoom: Mouse wheel<br/>
                        <Typography variant="h5">2-d Embedding</Typography>
                        Pan: Mouse left click<br/>
                        Zoom: Mouse wheel<br/>
                    </div>
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


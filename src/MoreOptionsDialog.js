import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import {debounce} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setChartOptions, setDialog} from './actions';


class MoreOptionsDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {labelFontSize: '', labelStrokeWidth: ''};
        this.onLabelFontSizeUpdate = debounce(this.onLabelFontSizeUpdate, 500);
        this.onLabelStrokeWidthUpdate = debounce(this.onLabelStrokeWidthUpdate, 500);
    }

    componentDidMount() {
        this.setState({
            labelFontSize: this.props.chartOptions.labelFontSize,
            labelStrokeWidth: this.props.chartOptions.labelStrokeWidth
        });
    }

    handleClose = () => {
        this.props.handleClose();
    };

    onLabelFontSizeUpdate = (value) => {
        if (!isNaN(value) && value > 0) {
            this.props.chartOptions.labelFontSize = value;
            this.props.handleChartOptions(this.props.chartOptions);
            this.setState({labelFontSize: value});
        }
    };

    onLabelFontSize = (event) => {
        this.setState({labelFontSize: event.target.value});
        this.onLabelFontSizeUpdate(event.target.value);
    };

    onLabelStrokeWidth = (event) => {
        this.setState({labelStrokeWidth: event.target.value});
        this.onLabelStrokeWidthUpdate(event.target.value);
    };

    onLabelStrokeWidthUpdate = (value) => {
        if (!isNaN(value) && value >= 0) {
            this.props.chartOptions.labelStrokeWidth = value;
            this.props.handleChartOptions(this.props.chartOptions);
            this.setState({labelStrokeWidth: value});
        }

    };

    render() {
        return (
            <Dialog
                open={true}
                onClose={this.handleClose}
                aria-labelledby="more-options-dialog-title"
                fullWidth={true}
                maxWidth={'sm'}
            >
                <DialogTitle id="more-options-dialog-title">More Options</DialogTitle>
                <DialogContent>

                    <TextField
                        value={this.state.labelFontSize}
                        onChange={this.onLabelFontSize}
                        margin="dense"
                        label="Label Font Size"
                        fullWidth
                    />

                    <TextField
                        value={this.state.labelStrokeWidth}
                        onChange={this.onLabelStrokeWidth}
                        margin="dense"
                        label="Label Shadow Size"
                        fullWidth
                    />

                </DialogContent>

            </Dialog>
        );
    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleChartOptions: value => {
            dispatch(setChartOptions(value));
        },
        handleClose: () => {
            dispatch(setDialog(null));
        },

    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(MoreOptionsDialog));


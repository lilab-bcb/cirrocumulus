import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';

import DialogTitle from '@material-ui/core/DialogTitle';
import TextField from '@material-ui/core/TextField';
import React from 'react';
import {connect} from 'react-redux';
import {setChartOptions, setDialog} from './actions';


class MoreOptionsDialog extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {labelFontSize: '', labelShadowSize: ''};
    }

    componentDidMount() {
        this.setState({
            labelFontSize: this.props.chartOptions.labelFontSize,
            labelShadowSize: this.props.chartOptions.labelShadowSize
        });
    }

    handleClose = () => {
        this.props.handleClose();
    };

    onLabelFontSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            let value = parseFloat(event.target.value);
            if (!isNaN(value) && value > 0) {
                this.props.chartOptions.labelFontSize = value;
                this.props.handleChartOptions(this.props.chartOptions);
                this.setState({labelFontSize: value});
            }
        }
    };

    onLabelFontSize = (event) => {
        this.setState({labelFontSize: event.target.value});
    };

    onLabelShadowSize = (event) => {
        this.setState({labelShadowSize: event.target.value});
    };

    onLabelShadowSizeKeyPress = (event) => {
        if (event.key === 'Enter') {
            let value = parseFloat(event.target.value);
            if (!isNaN(value) && value >= 0) {
                this.props.chartOptions.labelShadowSize = value;
                this.props.handleChartOptions(this.props.chartOptions);
                this.setState({labelShadowSize: value});
            }
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
                        onKeyPress={this.onLabelFontSizeKeyPress}
                        margin="dense"
                        label="Label Font Size"
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


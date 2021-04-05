import Divider from '@material-ui/core/Divider';
import React from 'react';

import {connect} from 'react-redux';
import {setPrimaryChartSize} from './actions';


class DraggableDivider extends React.PureComponent {


    onMouseDown = (event) => {
        this.dragging = true;
        this.clientY = event.clientY;
        this.primaryChartHeight = this.props.primaryChartSize.height;
        document.body.style.cursor = 'ns-resize';
        window.addEventListener('mousemove', this.onMouseMove);
        window.addEventListener('mouseup', this.onMouseUp);
    };

    onMouseUp = (event) => {
        if (this.dragging) {
            window.removeEventListener('mousemove', this.onMouseMove);
            window.removeEventListener('mouseup', this.onMouseUp);
            document.body.style.cursor = null;
        }
        this.dragging = false;
    };

    onMouseMove = (event) => {
        if (this.dragging) {
            const primaryChartSize = this.props.primaryChartSize;
            const delta = this.clientY - event.clientY;
            this.props.handlePrimaryChartSize({
                width: primaryChartSize.width,
                height: Math.max(50, this.primaryChartHeight - delta)
            });
        }
    };

    render() {
        const {activeFeature} = this.props;

        return (
            <>
                {activeFeature && <div style={{
                    height: 10, cursor: 'ns-resize', display: 'flex',
                    alignItems: 'center', justifyContent: 'center'
                }}
                                       onMouseDown={this.onMouseDown}>
                    <Divider style={{width: '100%'}}/>
                </div>}
            </>
        );
    }
}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        primaryChartSize: state.primaryChartSize
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handlePrimaryChartSize: value => {
            dispatch(setPrimaryChartSize(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DraggableDivider));


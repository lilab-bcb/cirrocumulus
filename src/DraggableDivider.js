import Divider from '@mui/material/Divider';
import React, {useRef} from 'react';

import {connect} from 'react-redux';
import {setPrimaryChartSize} from './actions';


function DraggableDivider(props) {
    const dragging = useRef(false);
    const clientY = useRef(-1);
    const primaryChartHeight = useRef(-1);

    function onMouseDown(event) {
        dragging.current = true;
        clientY.current = event.clientY;
        primaryChartHeight.current = props.primaryChartSize.height;
        document.body.style.cursor = 'ns-resize';
        window.addEventListener('mousemove', onMouseMove);
        window.addEventListener('mouseup', onMouseUp);
    }

    function onMouseUp(event) {
        if (dragging.current) {
            window.removeEventListener('mousemove', onMouseMove);
            window.removeEventListener('mouseup', onMouseUp);
            document.body.style.cursor = null;
        }
        dragging.current = false;
    }

    function onMouseMove(event) {
        if (dragging.current) {
            const primaryChartSize = props.primaryChartSize;
            const delta = clientY.current - event.clientY;
            props.handlePrimaryChartSize({
                width: primaryChartSize.width,
                height: Math.max(50, primaryChartHeight.current - delta)
            });
        }
    }


    const {activeFeature} = props;

    return (
        <div style={{
            height: 10, cursor: 'ns-resize',
            display: activeFeature ? 'flex' : 'none',
            alignItems: 'center', justifyContent: 'center'
        }}
             onMouseDown={onMouseDown}>
            <Divider style={{width: '100%'}}/>
        </div>
    );

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


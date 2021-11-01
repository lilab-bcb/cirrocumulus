import Divider from '@mui/material/Divider';
import React, {useRef} from 'react';

import {connect} from 'react-redux';
import {setDragDivider} from './actions';


function DraggableDivider(props) {
    const dragging = useRef(false);
    const clientY = useRef(-1);
    const primaryChartSizeHeightMouseDown = useRef();
    const {activeFeature, onDragDivider, primaryChartSizeHeight} = props;

    function onMouseDown(event) {
        dragging.current = true;
        clientY.current = event.clientY;
        primaryChartSizeHeightMouseDown.current = primaryChartSizeHeight;
        document.body.style.cursor = 'ns-resize';
        event.stopPropagation();
        window.addEventListener('mousemove', onMouseMove);
        window.addEventListener('mouseup', onMouseUp);
    }

    function onMouseUp(event) {
        if (dragging.current) {
            window.removeEventListener('mousemove', onMouseMove);
            window.removeEventListener('mouseup', onMouseUp);
            document.body.style.cursor = null;
            event.stopPropagation();
        }
        dragging.current = false;
    }

    function onMouseMove(event) {
        if (dragging.current) {
            const delta = clientY.current - event.clientY;
            const height = Math.max(50, primaryChartSizeHeightMouseDown.current - delta);
            onDragDivider(height);
            event.stopPropagation();
        }
    }


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
        primaryChartSizeHeight: state.panel.primaryChartSize.height
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onDragDivider: value => {
            dispatch(setDragDivider(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(DraggableDivider));


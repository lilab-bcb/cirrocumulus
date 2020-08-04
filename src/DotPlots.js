import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotSortOrder} from './actions';
import DotPlotCanvas from './DotPlotCanvas';


class DotPlots extends React.PureComponent {
    render() {
        const dotPlotData = this.props.dotPlotData;
        const selectedDotPlotData = this.props.selectedDotPlotData;
        if (dotPlotData.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        return <div>
            <div>{dotPlotData.map((data) => {
                return <DotPlotCanvas onSortOrderChanged={this.props.onSortOrderChanged} key={data.name} data={data}/>;
            })}</div>
            <div>{selectedDotPlotData.map((data) => {
                return <DotPlotCanvas onSortOrderChanged={this.props.onSortOrderChanged} key={data.name + '-selection'}
                                      data={data}/>;
            })}</div>
        </div>;
    }
}

const mapStateToProps = state => {
    return {
        dotPlotData: state.dotPlotData,
        selectedDotPlotData: state.selectedDotPlotData
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onSortOrderChanged: (payload) => {
            dispatch(setDotPlotSortOrder(payload));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DotPlots));


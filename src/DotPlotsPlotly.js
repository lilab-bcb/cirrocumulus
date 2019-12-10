import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotSortOrder} from './actions';
import DotPlot from './DotPlot';


class DotPlotsPlotly extends React.PureComponent {
    render() {
        return (<div>{this.props.dotPlotData.filter(data => data.active).map((data, i) => {
            return <DotPlot onSortOrderChanged={this.props.onSortOrderChanged} key={data.name} data={data}/>;
        })}</div>);
    }
}

const mapStateToProps = state => {
    return {
        dotPlotData: state.dotPlotData
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
)(DotPlotsPlotly));


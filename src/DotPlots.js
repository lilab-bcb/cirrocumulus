import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotSortOrder} from './actions';
import DotPlotCanvas from './DotPlotCanvas';


class DotPlots extends React.PureComponent {
    render() {
        const activeDotPlots = this.props.dotPlotData.filter(data => data.active);
        if (activeDotPlots.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        return (<div>{activeDotPlots.map((data, i) => {
            return <DotPlotCanvas onSortOrderChanged={this.props.onSortOrderChanged} key={data.name} data={data}/>;
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
)(DotPlots));


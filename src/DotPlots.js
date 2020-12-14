import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotInterpolator, setDotPlotOptions} from './actions';
import {DotPlotGroup} from './DotPlotGroup';
import {splitSearchTokens} from './util';

class DotPlots extends React.PureComponent {
    render() {
        const {
            chartOptions,
            dotPlotData,
            dotPlotInterpolator,
            categoricalNames,
            dotPlotOptions,
            handleInterpolator,
            onDotPlotOptions,
            selectedDotPlotData,
            searchTokens
        } = this.props;

        if (dotPlotData.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        const textColor = chartOptions.darkMode ? 'white' : 'black';
        const splitTokens = splitSearchTokens(searchTokens);
        const dimension = splitTokens.obsCat.join('-');
        let renamedCategories = categoricalNames[dimension] || {}; // TODO rename multiple
        return <DotPlotGroup
            dotPlotData={dotPlotData}
            selectedData={selectedDotPlotData}
            interpolator={dotPlotInterpolator}
            handleInterpolator={handleInterpolator}
            onDotPlotOptions={onDotPlotOptions}
            dotPlotOptions={dotPlotOptions}
            renamedCategories={renamedCategories}
            textColor={textColor}/>;

    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        dotPlotInterpolator: state.dotPlotInterpolator,
        dotPlotData: state.dotPlotData,
        dotPlotOptions: state.dotPlotOptions,
        selectedDotPlotData: state.selectedDotPlotData,
        categoricalNames: state.categoricalNames,
        searchTokens: state.searchTokens
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onDotPlotOptions: (payload) => {
            dispatch(setDotPlotOptions(payload));
        },
        handleInterpolator: value => {
            dispatch(setDotPlotInterpolator(value));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DotPlots));


import {groupBy} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setDistributionPlotInterpolator, setDistributionPlotOptions} from './actions';
import DistributionGroup from './DistributionGroup';


class DistributionPlots extends React.PureComponent {
    render() {
        const {
            cachedData,
            categoricalNames,
            chartOptions,
            distributionData,
            distributionPlotOptions,
            dotPlotInterpolator,
            globalFeatureSummary,
            handleInterpolator,
            onDistributionPlotOptions,
            selectedDistributionData
        } = this.props;

        if (distributionData.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        const textColor = chartOptions.darkMode ? 'white' : 'black';
        let dimension2data = groupBy(distributionData, 'dimension');
        let dimension2selecteddata = groupBy(selectedDistributionData, 'dimension');

        return <React.Fragment>{Object.keys(dimension2data).map(dimension => {
            const data = dimension2data[dimension];

            return <DistributionGroup key={dimension}
                                      cachedData={cachedData}
                                      distributionData={data}
                                      globalFeatureSummary={globalFeatureSummary}
                                      selectedData={dimension2selecteddata[dimension]}
                                      interpolator={dotPlotInterpolator}
                                      handleInterpolator={handleInterpolator}
                                      onDistributionPlotOptions={onDistributionPlotOptions}
                                      distributionPlotOptions={distributionPlotOptions}
                                      categoricalNames={categoricalNames}
                                      textColor={textColor}/>;
        })}</React.Fragment>;
    }


}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        distributionData: state.distributionData,
        distributionPlotOptions: state.distributionPlotOptions,
        dotPlotInterpolator: state.dotPlotInterpolator,
        globalFeatureSummary: state.globalFeatureSummary,
        selectedDistributionData: state.selectedDistributionData
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onDistributionPlotOptions: (payload) => {
            dispatch(setDistributionPlotOptions(payload));
        },
        handleInterpolator: value => {
            dispatch(setDistributionPlotInterpolator(value));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DistributionPlots));


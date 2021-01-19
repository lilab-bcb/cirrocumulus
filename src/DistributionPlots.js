import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
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
            embeddingData,
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
            const categoryColorScales = [];
            data[0].dimensions.forEach(dimension => {
                let found = false;
                for (let i = 0; i < embeddingData.length; i++) {
                    if (dimension === embeddingData[i].name) {
                        categoryColorScales.push(embeddingData[i].colorScale); // TODO make color scale independent of embedding
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    categoryColorScales.push(scaleOrdinal(schemeCategory10)); // TODO make color scale independent of embedding
                }

            });
            return <DistributionGroup key={dimension}
                                      cachedData={cachedData}
                                      categoryColorScales={categoryColorScales}
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
        embeddingData: state.embeddingData,
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


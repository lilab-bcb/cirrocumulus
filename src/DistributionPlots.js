import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import {groupBy} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setDistributionPlotInterpolator, setDistributionPlotOptions} from './actions';
import DistributionGroup from './DistributionGroup';


class DistributionPlots extends React.PureComponent {

    onInterpolator = (value) => {
        const scale = this.props.distributionPlotInterpolator.scale;
        value.scale = scale;
        this.props.onInterpolator(value);
    };

    onColorScalingChange = (value) => {
        const distributionPlotInterpolator = this.props.distributionPlotInterpolator;
        distributionPlotInterpolator.scale = value;
        this.props.onInterpolator(Object.assign({}, distributionPlotInterpolator));
    };

    render() {
        const {
            cachedData,
            categoricalNames,
            chartOptions,
            dataset,
            distributionData,
            distributionPlotOptions,
            distributionPlotInterpolator,
            embeddingData,
            globalFeatureSummary,
            onDistributionPlotOptions,
            selectedDistributionData,
            setTooltip
        } = this.props;

        if (distributionData.length === 0) {
            return null;
        }
        const textColor = chartOptions.darkMode ? 'white' : 'black';
        let dimension2data = groupBy(distributionData, 'dimension');
        let dimension2selecteddata = groupBy(selectedDistributionData, 'dimension');

        return <>{Object.keys(dimension2data).map(dimension => {
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
                                      setTooltip={setTooltip}
                                      categoryColorScales={categoryColorScales}
                                      dataset={dataset}
                                      distributionData={data}
                                      globalFeatureSummary={globalFeatureSummary}
                                      selectedData={dimension2selecteddata[dimension]}
                                      interpolator={distributionPlotInterpolator}
                                      distributionPlotOptions={distributionPlotOptions}
                                      categoricalNames={categoricalNames}
                                      textColor={textColor}
                                      handleInterpolator={this.onInterpolator}
                                      onColorScalingChange={this.onColorScalingChange}
                                      onDistributionPlotOptions={onDistributionPlotOptions}/>;
        })}</>;
    }


}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        dataset: state.dataset,
        distributionData: state.distributionData,
        distributionPlotOptions: state.distributionPlotOptions,
        distributionPlotInterpolator: state.distributionPlotInterpolator,
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
        onInterpolator: value => {
            dispatch(setDistributionPlotInterpolator(value));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DistributionPlots));


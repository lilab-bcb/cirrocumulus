import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import {groupBy, partial} from 'lodash';
import React from 'react';
import {connect} from 'react-redux';
import {setDistributionPlotInterpolator, setDistributionPlotOptions} from './actions';
import DistributionGroup from './DistributionGroup';
import Typography from '@mui/material/Typography';
import {DIST_PLOT_OPTIONS, DISTRIBUTION_PLOT_INTERPOLATOR_OBJ} from './reducers';


function DistributionPlots(props) {


    function onInterpolator(dataType, value) {
        const existingInterpolator = props.distributionPlotInterpolator[dataType];
        value.scale = existingInterpolator.scale;
        props.distributionPlotInterpolator[dataType] = value;
        props.onInterpolator(Object.assign({}, props.distributionPlotInterpolator));
    }

    function onColorScalingChange(dataType, value) {
        const existingInterpolator = props.distributionPlotInterpolator[dataType];
        existingInterpolator.scale = value;
        props.onInterpolator(Object.assign({}, props.distributionPlotInterpolator));
    }

    function onDistributionPlotOptionsDataType(dataType, value) {
        props.distributionPlotOptions[dataType] = Object.assign({}, props.distributionPlotOptions[dataType], value);
        onDistributionPlotOptions(Object.assign({}, props.distributionPlotOptions));
    }

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
        selectedDistributionData
    } = props;
    const textColor = chartOptions.darkMode ? 'white' : 'black';
    const keys = Object.keys(distributionData);
    const typeToName = {
        X: 'Features',
        modules: 'Modules',
        obs: 'Observations'
    };
    keys.sort((a, b) => {
        return a.toLowerCase() - b.toLowerCase();
    });
    const xIndex = keys.indexOf('X'); // put features 1st
    if (xIndex !== -1) {
        keys.splice(xIndex, 1);
        keys.unshift('X');
    }
    return keys.map(key => {
        const name = typeToName[key] || key;
        if (distributionPlotOptions[key] == null) {
            distributionPlotOptions[key] = Object.assign({}, DIST_PLOT_OPTIONS, {chartType: 'heatmap'});
            distributionPlotInterpolator[key] = Object.assign({}, DISTRIBUTION_PLOT_INTERPOLATOR_OBJ);
        }
        let dimension2data = groupBy(distributionData[key], 'dimension');
        let dimension2selecteddata = groupBy(selectedDistributionData[key], 'dimension');
        return <div key={key}><Typography color="textPrimary"
                                          variant={"h5"}>{name}</Typography>{Object.keys(dimension2data).map(dimension => {
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
                                      showDotPlotOption={key === 'X'}
                                      categoryColorScales={categoryColorScales}
                                      dataset={dataset}
                                      distributionData={data}
                                      globalFeatureSummary={globalFeatureSummary}
                                      selectedData={dimension2selecteddata[dimension]}
                                      interpolator={distributionPlotInterpolator[key]}
                                      distributionPlotOptions={distributionPlotOptions[key]}
                                      categoricalNames={categoricalNames}
                                      textColor={textColor}
                                      handleInterpolator={partial(onInterpolator, key)}
                                      onColorScalingChange={partial(onColorScalingChange, key)}
                                      onDistributionPlotOptions={partial(onDistributionPlotOptionsDataType, key)}/>;
        })}</div>;
    });


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
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(DistributionPlots));


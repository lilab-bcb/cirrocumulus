import React from 'react';

import {connect} from 'react-redux';
import {
    getEmbeddingKey,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated,
    handleSelectedPoints
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import createPlotlyComponent from './factory';
import PlotUtil from './PlotUtil';

const Plot = createPlotlyComponent(window.Plotly);

class EmbeddingChartPlotly extends React.PureComponent {


    getPlots() {
        const activeTraces = this.props.data.filter(traceInfo => traceInfo.active);
        const {config, nObs, nObsSelected, embeddingChartSize, onDeselect, onSelect, globalFeatureSummary, featureSummary, datasetFilter, handleDimensionFilterUpdated, handleMeasureFilterUpdated} = this.props;
        let size = PlotUtil.getEmbeddingChartSize(activeTraces.length === 1 ? 1 : embeddingChartSize);
        return activeTraces.map(traceInfo => {
            if (size !== traceInfo.layout.width) {
                traceInfo.layout = Object.assign({}, traceInfo.layout);
                traceInfo.layout.width = size;
                traceInfo.layout.height = size;
            }
            if (traceInfo.data[0].marker.size === 0) {
                traceInfo.data[0].marker.size = 1e-8;
            }
            if (traceInfo.data[0].unselected.marker.size === 0) {
                traceInfo.data[0].unselected.marker.size = 1e-8;
            }
            const key = traceInfo.name + '_' + getEmbeddingKey(traceInfo.data[0].embedding);
            return (
                <div style={{display: 'inline-block', border: '1px solid LightGrey'}} key={key}><Plot
                    data={traceInfo.data}
                    layout={traceInfo.layout}
                    config={config}
                    onDeselect={onDeselect}
                    onSelected={onSelect}
                />
                    {traceInfo.continuous ?
                        <ColorSchemeLegendWrapper
                            width={186}
                            label={true}
                            height={40}
                            handleUpdate={handleMeasureFilterUpdated}
                            datasetFilter={datasetFilter}
                            scale={traceInfo.colorScale}
                            featureSummary={featureSummary}
                            globalFeatureSummary={globalFeatureSummary}
                            nObs={nObs}
                            nObsSelected={nObsSelected}
                            maxHeight={traceInfo.layout.height}
                            name={traceInfo.name}
                        /> :
                        <CategoricalLegend datasetFilter={datasetFilter}
                                           handleClick={handleDimensionFilterUpdated}
                                           name={traceInfo.name}
                                           scale={traceInfo.colorScale}
                                           maxHeight={traceInfo.layout.height - 24}
                                           clickEnabled={true}
                                           nObs={nObs}
                                           nObsSelected={nObsSelected}
                                           globalFeatureSummary={globalFeatureSummary}
                                           featureSummary={featureSummary}/>}</div>);
        });
    }

    render() {
        return this.getPlots();
    }
}

const mapStateToProps = state => {
    return {
        data: state.embeddingData,
        numberOfBins: state.numberOfBins,
        binValues: state.binValues,
        embeddingChartSize: state.embeddingChartSize,
        config: state.plotConfig,
        datasetFilter: state.datasetFilter,
        featureSummary: state.featureSummary,
        nObs: state.dataset.nObs,
        nObsSelected: state.selection.count,
        globalFeatureSummary: state.globalFeatureSummary
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleDimensionFilterUpdated: (e) => {
            dispatch(handleDimensionFilterUpdated(e));
        },
        handleMeasureFilterUpdated: (e) => {
            dispatch(handleMeasureFilterUpdated(e));
        },
        onSelect: (e) => {
            dispatch(handleSelectedPoints(e));
        },
        onDeselect: () => {
            dispatch(handleSelectedPoints(null));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingChartPlotly));


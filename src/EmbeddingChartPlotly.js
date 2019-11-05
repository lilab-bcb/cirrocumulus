import React from 'react';

import {connect} from 'react-redux';
import {handleLegendClick, handleSelectedPoints} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import createPlotlyComponent from './factory';
import PlotUtil from './PlotUtil';

const Plot = createPlotlyComponent(window.Plotly);

class EmbeddingChartPlotly extends React.PureComponent {

    constructor(props) {
        super(props);
    }

    getKey(traceInfo) {
        let key = traceInfo.name;
        if (this.props.binValues) {
            key += '-' + this.props.numberOfBins;
        }
        return key;
    }

    getPlots() {
        const activeTraces = this.props.data.filter(traceInfo => traceInfo.active);
        const {config, embeddingChartSize, onDeselect, onSelect, featureSummary, nObs, categoricalFilter, handleLegendClick} = this.props;
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
            return (
                <div style={{display: 'inline-block', border: '1px solid LightGrey'}} key={traceInfo.name}><Plot
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
                            scale={traceInfo.colorScale}
                            featureSummary={featureSummary}
                            maxHeight={traceInfo.layout.height}
                            name={traceInfo.name}
                        /> :
                        <CategoricalLegend categoricalFilter={categoricalFilter}
                                           handleClick={handleLegendClick}
                                           name={traceInfo.name}
                                           scale={traceInfo.colorScale}
                                           maxHeight={traceInfo.layout.height - 24}
                                           clickEnabled={true}
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
        categoricalFilter: state.categoricalFilter,
        featureSummary: state.featureSummary
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleLegendClick: (e) => {
            dispatch(handleLegendClick(e));
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


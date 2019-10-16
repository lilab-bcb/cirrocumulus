import React from 'react';
import createPlotlyComponent from 'react-plotly.js/factory';

import {connect} from 'react-redux';
import {handleLegendClick, handleSelectedPoints} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';

const Plot = createPlotlyComponent(window.Plotly);

class EmbeddingChartPlotly extends React.PureComponent {


    getPlots() {
        return this.props.data.filter(traceInfo => traceInfo.active).map(traceInfo => {
            return (<div style={{display: 'inline-block', border: '1px solid LightGrey'}} key={traceInfo.name}><Plot
                key={traceInfo.name}
                data={traceInfo.data}
                layout={traceInfo.layout}
                config={this.props.config}
                onDeselect={this.props.onDeselect}
                onSelected={this.props.onSelect}
            />
                {traceInfo.continuous ?
                    <ColorSchemeLegendWrapper style={{display: 'block', marginLeft: 'auto', marginRight: 'auto'}}
                                              width={300}
                                              label={true} height={40}
                                              scale={traceInfo.colorScale}
                                              name={traceInfo.name}
                                              selectedValueCounts={this.props.selectedValueCounts}/> :
                    <CategoricalLegend legendVisibility={this.props.legendVisibility}
                                       handleClick={this.props.handleLegendClick} name={traceInfo.name}
                                       scale={traceInfo.colorScale}
                                       selectedValueCounts={this.props.selectedValueCounts}/>}</div>);
        });
    }

    render() {
        return this.getPlots();
    }
}

const mapStateToProps = state => {
    return {
        data: state.embeddingData,
        config: state.plotConfig,
        legendVisibility: state.legendVisibility,
        selectedValueCounts: state.selectedValueCounts
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
        },
        onZoom: () => {

        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingChartPlotly));


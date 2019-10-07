import {scaleSequential} from 'd3-scale';
import React from 'react';

import {connect} from 'react-redux';
import ColorSchemeLegend from './ColorSchemeLegend';
import PlotUtil from './PlotUtil';
import createPlotlyComponent from 'react-plotly.js/factory';

const Plot = createPlotlyComponent(window.Plotly);

class EmbeddingChartPlotly extends React.PureComponent {

    onRelayout = (eventdata) => {
        console.log('onRelayout')
        if (eventdata.xaxis) {
            if (eventdata.xaxis.autorange) {

            }
        }
        // let xaxisRange = eventdata.xaxis.range;
        // let yaxisRange = eventdata.yaxis.range;
        // let zaxisRange = eventdata.zaxis ? eventdata.zaxis.range : null;
        // console.log(xaxisRange, yaxisRange);

    };

    fixLegend = (figure, graphDiv) => {
        PlotUtil.fixLegend(graphDiv);
    };

    getScale(traceInfo) {
        if (traceInfo.data[0].domain != null) {
            return scaleSequential(this.props.interpolator.value).domain(traceInfo.data[0].domain);
        }
    }

    getPlots() {
        return this.props.data.filter(traceInfo => traceInfo.active).map(traceInfo => {
            return <div style={{display: 'inline-block', border: '1px solid LightGrey'}} key={traceInfo.name}><Plot
                key={traceInfo.name}
                data={traceInfo.data}
                layout={traceInfo.layout}
                config={this.props.config}
                onRelayout={this.onRelayout}
                onDeselect={this.props.onDeselect}
                onSelected={this.props.onSelect}
                onUpdate={this.fixLegend}
                onInitialized={this.fixLegend}
            /><ColorSchemeLegend style={{display: 'block', marginLeft: 'auto', marginRight: 'auto'}}
                                 width={300}
                                 label={true} height={40}
                                 scale={this.getScale(traceInfo)}/></div>;
        });
    }

    render() {
        return this.getPlots()
    }
}

const mapStateToProps = state => {
    return {
        data: state.embeddingData,
        config: state.plotConfig,
        interpolator: state.interpolator,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onSelect: (e) => {

        },
        onDeselect: () => {

        },
        onZoom: () => {

        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingChartPlotly));


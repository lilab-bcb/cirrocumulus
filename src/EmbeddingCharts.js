import React from 'react';

import {connect} from 'react-redux';
import {getEmbeddingKey} from './actions';
import EmbeddingChart from './EmbeddingChart';
import PlotUtil from './PlotUtil';


class EmbeddingCharts extends React.PureComponent {
    getPlots() {
        const embeddingChartSize = this.props.embeddingChartSize;
        const activeTraces = this.props.embeddingData.filter(traceInfo => traceInfo.active);

        return activeTraces.map(traceInfo => {
            const key = traceInfo.name + '_' + getEmbeddingKey(traceInfo.data[0].embedding);
            let size = PlotUtil.getEmbeddingChartSize(activeTraces.length === 1 ? 1 : embeddingChartSize);

            if (size !== traceInfo.layout.width) {
                traceInfo.layout = Object.assign({}, traceInfo.layout);
                traceInfo.layout.width = size;
                traceInfo.layout.height = size;
            }

            if (traceInfo.data[0].marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].marker.size = 1e-8;
            }
            if (traceInfo.data[0].unselected.marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].unselected.marker.size = 1e-8;
            }
            return <EmbeddingChart traceInfo={traceInfo} key={key}/>;
        });
    }

    render() {
        return this.getPlots();
    }
}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        embeddingChartSize: state.embeddingChartSize
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


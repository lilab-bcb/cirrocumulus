import React from 'react';

import {connect} from 'react-redux';
import {getEmbeddingKey, getTraceKey} from './actions';
import EmbeddingChart from './EmbeddingChart';
import {getChartSize} from './util';

const emptySet = new Set();

class EmbeddingCharts extends React.PureComponent {

    render() {
        const {primaryTraceKey, embeddingData, markerOpacity, unselectedMarkerOpacity, selection} = this.props;
        let primaryTraces = embeddingData.filter(traceInfo => getTraceKey(traceInfo) === primaryTraceKey);
        const primaryTrace = primaryTraces.length === 1 ? primaryTraces[0] : null;
        let userPoints = emptySet;
        if (primaryTrace) {
            const embedding = primaryTrace.embedding;
            const fullName = getEmbeddingKey(embedding);
            const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
            userPoints = chartSelection ? chartSelection.userPoints : emptySet;
        }

        if (primaryTrace == null) {
            const chartSize = getChartSize();
            return <div style={{height: chartSize.height}}></div>;
        }
        return (<EmbeddingChart
                markerOpacity={markerOpacity}
                unselectedMarkerOpacity={unselectedMarkerOpacity}
                traceInfo={primaryTrace}
                selection={userPoints}
                color={primaryTrace.colors}
                onGallery={this.props.onGallery}
            />
        );
    }

}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        markerOpacity: state.markerOpacity,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        selection: state.selection,
        primaryTraceKey: state.primaryTraceKey
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


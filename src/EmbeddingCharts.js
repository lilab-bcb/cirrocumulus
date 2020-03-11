import React from 'react';

import {connect} from 'react-redux';
import {getTraceKey, setPrimaryTraceKey} from './actions';
import EmbeddingChart from './EmbeddingChart';
import GalleryImage from './GalleryImage';
import PlotUtil from './PlotUtil';


class EmbeddingCharts extends React.PureComponent {


    onChartSelected = (traceInfo) => {
        this.props.handlePrimaryTraceKey(getTraceKey(traceInfo));
    };

    getPlots() {
        const {primaryTraceKey, embeddingData} = this.props;
        //  const embeddingChartSize = this.props.embeddingChartSize;
        const activeTraces = embeddingData.filter(traceInfo => traceInfo.active);
        let primaryTraces = embeddingData.filter(traceInfo => getTraceKey(traceInfo) === primaryTraceKey);
        if (primaryTraces.length === 0 && activeTraces.length > 0) {
            primaryTraces = [activeTraces[0]];
        }
        console.log('getPlots', activeTraces.length, primaryTraces.length);
        const singleChartSize = PlotUtil.getSingleEmbeddingChartSize();
        // let itemSize = Math.floor(singleChartSize / embeddingChartSize);
        activeTraces.forEach(traceInfo => {

            // if (itemSize !== traceInfo.layout.width) {
            //     traceInfo.layout = Object.assign({}, traceInfo.layout); // force re-render
            //     traceInfo.layout.width = itemSize;
            //     traceInfo.layout.height = itemSize;
            // }

            if (traceInfo.data[0].marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].marker.size = 1e-8;
            }
            if (traceInfo.data[0].unselected.marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].unselected.marker.size = 1e-8;
            }
        });

        if (primaryTraces.length > 0) {
            let traceInfo = primaryTraces[0];
            if (singleChartSize !== traceInfo.layout.width) {
                traceInfo.layout = Object.assign({}, traceInfo.layout); // force re-render
                traceInfo.layout.width = singleChartSize;
                traceInfo.layout.height = singleChartSize;

            }
        }

        return (<div>
            {primaryTraces.map(traceInfo => <EmbeddingChart
                style={{display: 'block', border: '1px solid LightGrey'}} traceInfo={traceInfo}
                key={getTraceKey(traceInfo)}/>)}
            {activeTraces.map(traceInfo => <GalleryImage
                style={{display: 'inline-block', border: '1px solid LightGrey'}}
                traceInfo={traceInfo}
                selected={primaryTraceKey === getTraceKey(traceInfo)}
                onSelect={this.onChartSelected}
                key={getTraceKey(traceInfo)}/>)}
        </div>);


    }

    render() {
        return this.getPlots();
    }
}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        embeddingChartSize: state.embeddingChartSize,
        primaryTraceKey: state.primaryTraceKey
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handlePrimaryTraceKey: (value) => {
            dispatch(setPrimaryTraceKey(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


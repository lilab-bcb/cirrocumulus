import React from 'react';

import {connect} from 'react-redux';
import {ScatterGL} from 'scatter-gl';
import {getEmbeddingKey, getTraceKey, setPrimaryTraceKey} from './actions';
import EmbeddingChart from './EmbeddingChart';
import GalleryImage from './GalleryImage';


class EmbeddingCharts extends React.PureComponent {


    constructor(props) {
        super(props);
        const containerElement = document.createElement('div');
        const ratio = window.devicePixelRatio;
        const size = Math.round(400 / ratio);
        containerElement.style.position = 'absolute';
        containerElement.style.left = '-9999px';
        containerElement.style.width = size + 'px';
        containerElement.style.height = size + 'px';
        document.body.appendChild(containerElement);
        this.scatterGL = new ScatterGL(containerElement, {
            renderMode: 'POINT',
            interactive: false,
            rotateOnStart: false,
            showLabelsOnHover: false
        });
        this.containerElement = containerElement;
        this.emptySet = new Set();
    }

    onChartSelected = (traceInfo) => {
        this.props.handlePrimaryTraceKey(getTraceKey(traceInfo));
    };

    getPlots() {
        const {primaryTraceKey, embeddingData, markerOpacity, unselectedMarkerOpacity, selection} = this.props;
        let activeTraces = embeddingData.filter(traceInfo => traceInfo.active);
        let primaryTraces = embeddingData.filter(traceInfo => getTraceKey(traceInfo) === primaryTraceKey);
        if (primaryTraces.length === 0 && activeTraces.length > 0) {
            primaryTraces = [activeTraces[0]];
        }
        activeTraces = activeTraces.filter(activeTrace => activeTrace.name !== '__count');

        const primaryTrace = primaryTraces[0];
        let userPoints = this.emptySet;
        if (primaryTrace) {
            const embedding = primaryTrace.embedding;
            const fullName = getEmbeddingKey(embedding);
            const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
            userPoints = chartSelection ? chartSelection.userPoints : this.emptySet;
        }
        return (<React.Fragment>
            {primaryTrace && <EmbeddingChart
                markerOpacity={markerOpacity}
                unselectedMarkerOpacity={unselectedMarkerOpacity}
                traceInfo={primaryTrace}
                selection={userPoints}
                color={primaryTrace.marker.color}
            />}
            {activeTraces.map(traceInfo => <GalleryImage
                traceInfo={traceInfo}
                color={traceInfo.marker.color}
                scatterGL={this.scatterGL}
                markerOpacity={markerOpacity}
                containerElement={this.containerElement}
                onSelect={this.onChartSelected}
                key={getTraceKey(traceInfo)}/>)}
        </React.Fragment>);


    }

    render() {
        return this.getPlots();
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
    return {
        handlePrimaryTraceKey: (value) => {
            dispatch(setPrimaryTraceKey(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


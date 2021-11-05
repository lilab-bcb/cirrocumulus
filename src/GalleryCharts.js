import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {getTraceKey, setActiveFeature, setEmbeddingData} from './actions';
import GalleryImage from './GalleryImage';
import {createScatterPlot} from './ThreeUtil';
import {findIndex} from 'lodash';
import {FEATURE_TYPE} from './util';

function createContainer(chartSize) {
    const containerElement = document.createElement('div');
    containerElement.style.position = 'absolute';
    containerElement.style.left = '-9999px';
    containerElement.style.width = chartSize + 'px';
    containerElement.style.height = chartSize + 'px';
    return containerElement;
}


class GalleryCharts extends React.PureComponent {
    constructor(props) {
        super(props);
        this.containerElement = createContainer(this.props.chartSize);
        document.body.appendChild(this.containerElement);
        this.scatterPlot = createScatterPlot(this.containerElement, window.ApplePaySession != null, false, false);
    }

    onChartSelected = (traceInfo) => {
        this.props.handleActiveFeature({
            name: traceInfo.name,
            type: traceInfo.featureType,
            embeddingKey: getTraceKey(traceInfo)
        });
        window.scrollTo(0, 0);
    };

    onSortEnd = (galleryTraces, e) => {
        const oldTrace = galleryTraces[e.oldIndex];
        const newTrace = galleryTraces[e.newIndex];
        const embeddingData = this.props.embeddingData;
        const oldIndex = findIndex(embeddingData, oldTrace);
        const newIndex = findIndex(embeddingData, newTrace);
        embeddingData.splice(oldIndex, 1);
        embeddingData.splice(newIndex, 0, oldTrace);
        this.props.handleEmbeddingData(embeddingData.slice());
    };


    render() {
        const {
            cachedData,
            categoricalNames,
            chartSize,
            chartOptions,
            embeddingData,
            embeddingLabels,
            markerOpacity,
            pointSize,
            primaryChartSize,
            searchTokens,
            selection,
            unselectedMarkerOpacity,
            unselectedPointSize
        } = this.props;

        if (this.containerElement.style.width !== this.props.chartSize + 'px') {
            document.body.removeChild(this.containerElement);
            this.containerElement = createContainer(this.props.chartSize);
            document.body.appendChild(this.containerElement);
            this.scatterPlot = createScatterPlot(this.containerElement, window.ApplePaySession, false, false);
        }

        const galleryTraces = embeddingData.filter(traceInfo => traceInfo.active);
        const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT && embeddingLabels.indexOf(item.id) !== -1).map(item => item.id);
        const SortableItem = sortableElement(({trace}) => <GalleryImage
            traceInfo={trace}
            obsCat={obsCat}
            cachedData={cachedData}
            scatterPlot={this.scatterPlot}
            markerOpacity={markerOpacity}
            chartOptions={chartOptions}
            pointSize={pointSize}
            unselectedPointSize={unselectedPointSize}
            primaryChartSize={primaryChartSize}
            chartSize={chartSize}
            categoricalNames={categoricalNames}
            unselectedMarkerOpacity={unselectedMarkerOpacity}
            selection={selection}
            containerElement={this.containerElement}
            onSelect={this.onChartSelected}
            key={getTraceKey(trace)}/>);

        const SortableList = sortableContainer(({items}) => {
            return (
                <ul style={{padding: 0, marginTop: 4, marginBottom: 0}}>
                    {items.map((trace, index) => (
                        <SortableItem key={getTraceKey(trace)} index={index} trace={trace}/>
                    ))}
                </ul>
            );
        });

        return (
            <SortableList
                distance={2}
                axis="xy" items={galleryTraces}
                onSortEnd={(e) => this.onSortEnd(galleryTraces, e)}/>
        );
    }
}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        chartSize: state.chartSize,
        embeddingData: state.embeddingData,
        embeddingLabels: state.embeddingLabels,
        markerOpacity: state.markerOpacity,
        pointSize: state.pointSize,
        primaryChartSize: state.panel.primaryChartSize,
        searchTokens: state.searchTokens,
        selection: state.selection,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        unselectedPointSize: state.unselectedPointSize
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleActiveFeature: (value) => {
            dispatch(setActiveFeature(value));
        },
        handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(GalleryCharts));


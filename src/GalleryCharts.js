import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {ScatterGL} from 'scatter-gl';
import {getTraceKey, setEmbeddingData, setPrimaryTraceKey} from './actions';
import GalleryImage from './GalleryImage';

class GalleryCharts extends React.PureComponent {


    constructor(props) {
        super(props);
        const containerElement = document.createElement('div');
        const size = 400;
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

    onSortEnd = (galleryTraces, e) => {
        galleryTraces[e.oldIndex].sortIndex = e.newIndex;
        galleryTraces[e.newIndex].sortIndex = e.oldIndex;
        this.props.handleEmbeddingData(this.props.embeddingData.slice());
    };

    render() {
        const {embeddingData, markerOpacity} = this.props;
        let galleryTraces = embeddingData.filter(traceInfo => traceInfo.active);

        for (let i = 0; i < galleryTraces.length; i++) {
            if (galleryTraces[i].sortIndex == null) {
                galleryTraces[i].sortIndex = i;
            }
        }
        galleryTraces.sort((a, b) => a.sortIndex - b.sortIndex);

        // const DragHandle = sortableHandle(() => <span>::</span>);

        const SortableItem = sortableElement(({traceInfo}) => <GalleryImage
            traceInfo={traceInfo}
            color={traceInfo.marker.color}
            scatterGL={this.scatterGL}
            markerOpacity={markerOpacity}
            containerElement={this.containerElement}
            onSelect={this.onChartSelected}
            key={getTraceKey(traceInfo)}/>);

        const SortableList = sortableContainer(({items}) => {
            return (
                <ul>
                    {items.map((traceInfo, index) => (
                        <SortableItem key={getTraceKey(traceInfo)} index={index} traceInfo={traceInfo}/>
                    ))}
                </ul>
            );
        });
        return (
            <SortableList distance={2}
                          axis="xy" items={galleryTraces}
                          onSortEnd={(e) => this.onSortEnd(galleryTraces, e)}/>);
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
        },
        handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(GalleryCharts));


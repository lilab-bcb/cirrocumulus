import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {ScatterGL} from 'scatter-gl';
import {getEmbeddingKey, getTraceKey, setEmbeddingData, setPrimaryTraceKey} from './actions';
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

    onSortEnd = (activeTraces, e) => {
        activeTraces[e.oldIndex].sortIndex = e.newIndex;
        activeTraces[e.newIndex].sortIndex = e.oldIndex;
        this.props.handleEmbeddingData(this.props.embeddingData.slice(0));
    };

    render() {
        const {primaryTraceKey, embeddingData, markerOpacity, unselectedMarkerOpacity, selection} = this.props;
        let activeTraces = embeddingData.filter(traceInfo => traceInfo.active);
        let primaryTraces = embeddingData.filter(traceInfo => getTraceKey(traceInfo) === primaryTraceKey);
        if (primaryTraces.length === 0 && activeTraces.length > 0) {
            primaryTraces = [activeTraces[0]];
        }
        activeTraces = activeTraces.filter(activeTrace => activeTrace.name !== '__count');
        for (let i = 0; i < activeTraces.length; i++) {
            if (activeTraces[i].sortIndex == null) {
                activeTraces[i].sortIndex = i;
            }
        }
        activeTraces.sort((a, b) => a.sortIndex - b.sortIndex);
        const primaryTrace = primaryTraces[0];
        let userPoints = this.emptySet;
        if (primaryTrace) {
            const embedding = primaryTrace.embedding;
            const fullName = getEmbeddingKey(embedding);
            const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
            userPoints = chartSelection ? chartSelection.userPoints : this.emptySet;
        }
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
        return (<React.Fragment>
            {primaryTrace && <EmbeddingChart
                markerOpacity={markerOpacity}
                unselectedMarkerOpacity={unselectedMarkerOpacity}
                traceInfo={primaryTrace}
                selection={userPoints}
                color={primaryTrace.marker.color}
            />}
            <SortableList distance={2}
                          axis="xy" items={activeTraces}
                          onSortEnd={(e) => this.onSortEnd(activeTraces, e)}/>;
        </React.Fragment>);
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
)(EmbeddingCharts));


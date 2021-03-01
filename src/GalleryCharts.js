import Divider from '@material-ui/core/Divider';
import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {getTraceKey, setActiveFeature, setEmbeddingData, setPrimaryChartSize} from './actions';
import GalleryImage from './GalleryImage';
import {createScatterPlot} from './ThreeUtil';
import {splitSearchTokens} from './util';

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
        galleryTraces[e.oldIndex].sortIndex = e.newIndex;
        galleryTraces[e.newIndex].sortIndex = e.oldIndex;
        this.props.handleEmbeddingData(this.props.embeddingData.slice());
    };


    onMouseDown = (event) => {
        this.dragging = true;
        this.clientY = event.clientY;
        this.primaryChartHeight = this.props.primaryChartSize.height;
        document.body.style.cursor = 'ns-resize';
        window.addEventListener('mousemove', this.onMouseMove);
        window.addEventListener('mouseup', this.onMouseUp);
    };

    onMouseUp = (event) => {
        if (this.dragging) {
            window.removeEventListener('mousemove', this.onMouseMove);
            window.removeEventListener('mouseup', this.onMouseUp);
            document.body.style.cursor = null;
        }
        this.dragging = false;
    };

    onMouseMove = (event) => {
        if (this.dragging) {
            const primaryChartSize = this.props.primaryChartSize;
            const delta = this.clientY - event.clientY;
            this.props.handlePrimaryChartSize({
                width: primaryChartSize.width,
                height: this.primaryChartHeight - delta
            });
        }
    };

    render() {
        const {
            activeFeature,
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
            unselectedMarkerOpacity
        } = this.props;
        if (this.containerElement.style.width !== this.props.chartSize + 'px') {
            document.body.removeChild(this.containerElement);
            this.containerElement = createContainer(this.props.chartSize);
            document.body.appendChild(this.containerElement);
            this.scatterPlot = createScatterPlot(this.containerElement, window.ApplePaySession, false, false);
        }
        const galleryTraces = embeddingData.filter(traceInfo => traceInfo.active);
        const obsCat = splitSearchTokens(searchTokens).obsCat.filter(item => embeddingLabels.indexOf(item) !== -1);
        // const DragHandle = sortableHandle(() => <span>::</span>);

        const SortableItem = sortableElement(({trace}) => <GalleryImage
            traceInfo={trace}
            obsCat={obsCat}
            cachedData={cachedData}
            scatterPlot={this.scatterPlot}
            markerOpacity={markerOpacity}
            chartOptions={chartOptions}
            pointSize={pointSize}
            chartSize={chartSize}
            categoricalNames={categoricalNames}
            primaryChartSize={primaryChartSize}
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
            <React.Fragment>
                {activeFeature && <div style={{
                    height: 10, cursor: 'ns-resize', display: 'flex',
                    alignItems: 'center', justifyContent: 'center'
                }}
                     onMouseDown={this.onMouseDown}>
                    <Divider style={{width: '100%'}}/>
                </div>}

                <SortableList
                    distance={2}
                    axis="xy" items={galleryTraces}
                    onSortEnd={(e) => this.onSortEnd(galleryTraces, e)}/></React.Fragment>
        );
    }
}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        chartSize: state.chartSize,
        embeddingData: state.embeddingData,
        embeddingLabels: state.embeddingLabels,
        markerOpacity: state.markerOpacity,
        pointSize: state.pointSize,
        primaryChartSize: state.primaryChartSize,
        searchTokens: state.searchTokens,
        selection: state.selection,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handlePrimaryChartSize: value => {
            dispatch(setPrimaryChartSize(value));
        },
        handleActiveFeature: (value) => {
            dispatch(setActiveFeature(value));
        },
        handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(GalleryCharts));


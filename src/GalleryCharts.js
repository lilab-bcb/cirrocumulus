import {Tooltip} from '@material-ui/core';
import Typography from '@material-ui/core/Typography';
import HelpOutlineIcon from '@material-ui/icons/HelpOutline';
import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {getTraceKey, setEmbeddingData, setPrimaryTraceKey} from './actions';
import GalleryImage from './GalleryImage';
import {createScatterPlot} from './ThreeUtil';

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
        this.scatterPlot = createScatterPlot(this.containerElement, window.ApplePaySession, false, false);
        this.emptySet = new Set();

    }

    onChartSelected = (traceInfo) => {
        this.props.handlePrimaryTraceKey(getTraceKey(traceInfo));
        window.scrollTo(0, 0);
    };

    onSortEnd = (galleryTraces, e) => {
        galleryTraces[e.oldIndex].sortIndex = e.newIndex;
        galleryTraces[e.newIndex].sortIndex = e.oldIndex;
        this.props.handleEmbeddingData(this.props.embeddingData.slice());
    };


    render() {
        let {categoricalNames, chartSize, embeddingData, primaryChartSize, markerOpacity, unselectedMarkerOpacity, pointSize, chartOptions, selection} = this.props;
        if (this.containerElement.style.width !== this.props.chartSize + 'px') {
            document.body.removeChild(this.containerElement);
            this.containerElement = createContainer(this.props.chartSize);
            document.body.appendChild(this.containerElement);
            this.scatterPlot = createScatterPlot(this.containerElement, window.ApplePaySession, false, false);
        }
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
            color={traceInfo.colors}
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
            key={getTraceKey(traceInfo)}/>);

        const SortableList = sortableContainer(({items}) => {
            return (
                <ul style={{padding: 0}}>
                    {items.map((traceInfo, index) => (
                        <SortableItem key={getTraceKey(traceInfo)} index={index} traceInfo={traceInfo}/>
                    ))}
                </ul>
            );
        });

        return (
            <React.Fragment><Typography variant="subtitle1">Gallery<Tooltip
                title="Drag gallery charts to reorder. Click chart to set primary view."><HelpOutlineIcon
                style={{verticalAlign: 'text-bottom'}}/></Tooltip></Typography>
                <SortableList
                    distance={2}
                    axis="xy" items={galleryTraces}
                    onSortEnd={(e) => this.onSortEnd(galleryTraces, e)}/></React.Fragment>);
    }
}

const mapStateToProps = state => {
    return {
        categoricalNames: state.categoricalNames,
        embeddingData: state.embeddingData,
        markerOpacity: state.markerOpacity,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        selection: state.selection,
        primaryTraceKey: state.primaryTraceKey,
        pointSize: state.pointSize,
        chartSize: state.chartSize,
        primaryChartSize: state.primaryChartSize,
        chartOptions: state.chartOptions
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


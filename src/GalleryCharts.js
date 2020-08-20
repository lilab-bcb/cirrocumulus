import {Tooltip} from '@material-ui/core';
import Typography from '@material-ui/core/Typography';
import HelpOutlineIcon from '@material-ui/icons/HelpOutline';
import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {getTraceKey, setEmbeddingData, setPrimaryTraceKey} from './actions';
import GalleryImage from './GalleryImage';
import {createScatterPlot} from './ThreeUtil';

class GalleryCharts extends React.PureComponent {
    constructor(props) {
        super(props);
        const containerElement = document.createElement('div');
        containerElement.style.position = 'absolute';
        containerElement.style.left = '-9999px';
        containerElement.style.width = props.chartSize + 'px';
        containerElement.style.height = props.chartSize + 'px';
        document.body.appendChild(containerElement);
        this.scatterPlot = createScatterPlot(containerElement, window.ApplePaySession);
        this.scatterPlot.interactive = false;
        this.containerElement = containerElement;
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

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.chartSize !== this.props.chartSize) {
            this.containerElement.style.width = this.props.chartSize + 'px';
            this.containerElement.style.height = this.props.chartSize + 'px';
            this.scatterPlot.resize(false);
        }
    }

    render() {
        const {chartSize, embeddingData, markerOpacity, unselectedMarkerOpacity, pointSize, chartOptions, selection} = this.props;
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
                title="Drag gallery charts to reorder. Click title to set primary view."><HelpOutlineIcon
                style={{verticalAlign: 'text-bottom'}}/></Tooltip></Typography>
                <SortableList
                    distance={2}
                    axis="xy" items={galleryTraces}
                    onSortEnd={(e) => this.onSortEnd(galleryTraces, e)}/></React.Fragment>);
    }
}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        markerOpacity: state.markerOpacity,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        selection: state.selection,
        primaryTraceKey: state.primaryTraceKey,
        pointSize: state.pointSize,
        chartSize: state.chartSize,
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


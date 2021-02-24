import {Tooltip} from '@material-ui/core';
import Typography from '@material-ui/core/Typography';
import HelpOutlineIcon from '@material-ui/icons/HelpOutline';
import React from 'react';

import {connect} from 'react-redux';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {getTraceKey, setActiveFeature, setEmbeddingData} from './actions';
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
        this.emptySet = new Set();

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


    render() {
        let {
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

        const SortableItem = sortableElement(({traceInfo}) => <GalleryImage
            traceInfo={traceInfo}
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
            <React.Fragment><Typography color={"textSecondary"} variant="subtitle1">Gallery<Tooltip
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


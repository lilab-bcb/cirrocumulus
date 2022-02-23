import React, {useEffect, useRef, useState} from 'react';

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


function GalleryCharts(props) {
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
        unselectedPointSize,
        handleActiveFeature,
        handleEmbeddingData
    } = props;

    const scatterPlotRef = useRef();
    const containerElementRef = useRef();
    const [forceUpdate, setForceUpdate] = useState(false);

    function onChartSelected(trace) {
        handleActiveFeature({
            name: trace.name, type: trace.featureType, embeddingKey: getTraceKey(trace)
        });
        window.scrollTo(0, 0);
    }

    function onSortEnd(galleryTraces, e) {
        const oldTrace = galleryTraces[e.oldIndex];
        const newTrace = galleryTraces[e.newIndex];
        const oldIndex = findIndex(embeddingData, oldTrace);
        const newIndex = findIndex(embeddingData, newTrace);
        embeddingData.splice(oldIndex, 1);
        embeddingData.splice(newIndex, 0, oldTrace);
        handleEmbeddingData(embeddingData.slice());
    }

    useEffect(() => {
        containerElementRef.current = createContainer(chartSize);
        document.body.appendChild(containerElementRef.current);
        scatterPlotRef.current = createScatterPlot(containerElementRef.current, window.ApplePaySession, false, false);
        const canvas = containerElementRef.current.querySelector('canvas');

        function webglcontextlost(e) {
            console.log('gallery - lost webgl context');
            e.preventDefault();
        }

        function webglcontextrestored(e) {
            console.log('gallery - restored webgl context');
            e.preventDefault();
            setForceUpdate(c => !c);
        }

        canvas.addEventListener('webglcontextlost', webglcontextlost);
        canvas.addEventListener('webglcontextrestored', webglcontextrestored);
        return () => {
            canvas.removeEventListener('webglcontextlost', webglcontextlost);
            canvas.removeEventListener('webglcontextrestored', webglcontextrestored);
            if (scatterPlotRef.current) {
                scatterPlotRef.current.dispose();
            }
            if (containerElementRef.current) {
                document.body.removeChild(containerElementRef.current);
            }
        };
    }, [chartSize, containerElementRef, scatterPlotRef]);


    const galleryTraces = embeddingData.filter(trace => trace.active);
    const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT && embeddingLabels.indexOf(item.id) !== -1).map(item => item.id);
    const SortableItem = sortableElement(({trace}) => <GalleryImage
        trace={trace}
        obsCat={obsCat}
        cachedData={cachedData}
        scatterPlot={scatterPlotRef.current}
        markerOpacity={markerOpacity}
        chartOptions={chartOptions}
        pointSize={pointSize}
        unselectedPointSize={unselectedPointSize}
        primaryChartSize={primaryChartSize}
        chartSize={chartSize}
        categoricalNames={categoricalNames}
        unselectedMarkerOpacity={unselectedMarkerOpacity}
        selection={selection}
        containerElement={containerElementRef.current}
        onSelect={onChartSelected}
        key={getTraceKey(trace)}/>);

    const SortableList = sortableContainer(({items}) => {
        return (<ul style={{padding: 0, marginTop: 4, marginBottom: 0}}>
            {items.map((trace, index) => (<SortableItem key={getTraceKey(trace)} index={index} trace={trace}/>))}
        </ul>);
    });

    return (<SortableList
        distance={2}
        axis="xy" items={galleryTraces}
        onSortEnd={(e) => onSortEnd(galleryTraces, e)}/>);
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
        }, handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        }
    };
};

export default (connect(mapStateToProps, mapDispatchToProps)(GalleryCharts));


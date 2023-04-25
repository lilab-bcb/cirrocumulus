import React, {memo, useEffect, useRef, useState} from 'react';
import {
  closestCenter,
  DndContext,
  PointerSensor,
  useSensor,
  useSensors,
} from '@dnd-kit/core';
import {
  arrayMove,
  rectSwappingStrategy,
  SortableContext,
  useSortable,
} from '@dnd-kit/sortable';
import {connect} from 'react-redux';

import {getTraceKey, setActiveFeature, setEmbeddingData} from './actions';
import GalleryImage from './GalleryImage';
import {createScatterPlot} from './ThreeUtil';
import {FEATURE_TYPE} from './util';
import {CSS} from '@dnd-kit/utilities';
import {find, findIndex} from 'lodash';

function SortableItem(props) {
  const {attributes, listeners, setNodeRef, transform, transition, isDragging} =
    useSortable({id: props.id});

  const style = {
    transform: CSS.Transform.toString(transform),
    transition,
    display: 'inline-block',
    zIndex: isDragging ? 1000 : 0,
    position: 'relative',
  };
  const {
    trace,
    obsCat,
    cachedData,
    displayName,
    markerOpacity,
    chartOptions,
    chartSize,
    pointSize,
    unselectedPointSize,
    primaryChartSize,
    categoricalNames,
    unselectedMarkerOpacity,
    selection,
    scatterPlot,
    containerElement,
    onSelect,
  } = props;

  return (
    <div ref={setNodeRef} style={style} {...attributes} {...listeners}>
      <MemoGalleryImage
        trace={trace}
        obsCat={obsCat}
        cachedData={cachedData}
        scatterPlot={scatterPlot}
        markerOpacity={markerOpacity}
        displayName={displayName}
        chartOptions={chartOptions}
        pointSize={pointSize}
        unselectedPointSize={unselectedPointSize}
        primaryChartSize={primaryChartSize}
        chartSize={chartSize}
        categoricalNames={categoricalNames}
        unselectedMarkerOpacity={unselectedMarkerOpacity}
        selection={selection}
        containerElement={containerElement}
        onSelect={onSelect}
        key={getTraceKey(trace)}
      />
    </div>
  );
}

const MemoGalleryImage = memo(GalleryImage);

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
    activeFeature,
    cachedData,
    categoricalNames,
    chartSize,
    chartOptions,
    embeddingData,
    embeddingLabels,
    jobResults,
    markerOpacity,
    pointSize,
    primaryChartSize,
    searchTokens,
    selection,
    tab,
    unselectedMarkerOpacity,
    unselectedPointSize,
    handleActiveFeature,
    handleEmbeddingData,
  } = props;

  const scatterPlotRef = useRef();
  const containerElementRef = useRef();
  const galleryTraces = embeddingData.filter((trace) => trace.active);

  function updateActiveTrace(e, index) {
    const trace = galleryTraces[index % galleryTraces.length];
    handleActiveFeature({
      name: trace.name,
      type: trace.featureType,
      embeddingKey: getTraceKey(trace),
    });
    window.scrollTo(0, 0);
    e.preventDefault();
    e.stopPropagation();
  }

  function handleKeyUp(e) {
    if (tab !== 'embedding') {
      return;
    }
    const tagName = e.target.tagName;
    const isInputField =
      tagName === 'INPUT' || tagName === 'SELECT' || tagName === 'TEXTAREA';
    if (!isInputField && galleryTraces.length > 0) {
      // 1-9
      if (e.key == '[' || e.key == ']') {
        // previous, next
        const previous = e.key == '[';
        let index = findIndex(
          galleryTraces,
          (trace) => getTraceKey(trace) === activeFeature.embeddingKey,
        );
        if (index < 0) {
          console.log('not found');
          return;
        }
        if (previous) {
          index--;
          if (index < 0) {
            index = galleryTraces.length - 1;
          }
        } else {
          index++;
          if (index === galleryTraces.length) {
            index = 0;
          }
        }
        updateActiveTrace(e, index);
      } else {
        const key = parseInt(e.key);
        if (key >= 1 && key <= 9) {
          updateActiveTrace(e, key - 1);
        }
      }
    }
  }

  const [forceUpdate, setForceUpdate] = useState(false);
  useEffect(() => {
    window.addEventListener('keypress', handleKeyUp);
    return () => {
      window.removeEventListener('keypress', handleKeyUp);
    };
  }, [handleKeyUp]);
  const sensors = useSensors(
    useSensor(PointerSensor, {
      activationConstraint: {
        distance: 1,
      },
    }),
  );

  function onChartSelected(trace) {
    handleActiveFeature({
      name: trace.name,
      type: trace.featureType,
      embeddingKey: getTraceKey(trace),
    });
    window.scrollTo(0, 0);
  }

  useEffect(() => {
    containerElementRef.current = createContainer(chartSize);
    document.body.appendChild(containerElementRef.current);
    scatterPlotRef.current = createScatterPlot(
      containerElementRef.current,
      window.ApplePaySession,
      false,
      false,
    );
    const canvas = containerElementRef.current.querySelector('canvas');

    function webglcontextlost(e) {
      console.log('gallery - lost webgl context');
      e.preventDefault();
    }

    function webglcontextrestored(e) {
      console.log('gallery - restored webgl context');
      e.preventDefault();
      setForceUpdate((c) => !c);
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

  const obsCat = searchTokens
    .filter(
      (item) =>
        item.type === FEATURE_TYPE.OBS_CAT &&
        embeddingLabels.indexOf(item.id) !== -1,
    )
    .map((item) => item.id);

  function dragEnd(event) {
    const {active, over} = event;
    if (active.id !== over.id) {
      const ids = embeddingData.map((trace) => getTraceKey(trace));
      const oldIndex = ids.indexOf(active.id);
      const newIndex = ids.indexOf(over.id);
      handleEmbeddingData(arrayMove(embeddingData, oldIndex, newIndex));
    }
  }

  function getDisplayName(trace) {
    let displayName;
    if (trace.featureType === FEATURE_TYPE.JOB_RESULT) {
      const val = find(jobResults, (jobResult) => jobResult.id === trace.name);
      displayName = val ? val.name : '';
    } else {
      displayName = trace.name === '__count' ? '' : trace.name;
    }
    return displayName;
  }

  return (
    <DndContext
      sensors={sensors}
      collisionDetection={closestCenter}
      onDragEnd={dragEnd}
    >
      <SortableContext
        strategy={rectSwappingStrategy}
        items={galleryTraces.map((trace) => getTraceKey(trace))}
      >
        {galleryTraces.map((trace) => (
          <SortableItem
            key={getTraceKey(trace)}
            id={getTraceKey(trace)}
            trace={trace}
            obsCat={obsCat}
            cachedData={cachedData}
            displayName={getDisplayName(trace)}
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
          />
        ))}
      </SortableContext>
    </DndContext>
  );
}

const mapStateToProps = (state) => {
  return {
    activeFeature: state.activeFeature,
    cachedData: state.cachedData,
    categoricalNames: state.categoricalNames,
    chartOptions: state.chartOptions,
    chartSize: state.chartSize,
    embeddingData: state.embeddingData,
    embeddingLabels: state.embeddingLabels,
    jobResults: state.jobResults,
    markerOpacity: state.markerOpacity,
    pointSize: state.pointSize,
    primaryChartSize: state.panel.primaryChartSize,
    searchTokens: state.searchTokens,
    selection: state.selection,
    tab: state.tab,
    unselectedMarkerOpacity: state.unselectedMarkerOpacity,
    unselectedPointSize: state.unselectedPointSize,
  };
};
const mapDispatchToProps = (dispatch) => {
  return {
    handleActiveFeature: (value) => {
      dispatch(setActiveFeature(value));
    },
    handleEmbeddingData: (value) => {
      dispatch(setEmbeddingData(value));
    },
  };
};

export default connect(mapStateToProps, mapDispatchToProps)(GalleryCharts);

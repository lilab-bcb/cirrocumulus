import {Tooltip} from '@mui/material';
import Box from '@mui/material/Box';
import CircularProgress from '@mui/material/CircularProgress';
import Typography from '@mui/material/Typography';
import React, {useEffect, useRef, useState} from 'react';
import {drawEmbeddingImage, getSpotRadius} from './ImageChart';
import {drawLabels, getVisualizer} from './ScatterChartThree';
import {
    getCategoryLabelsPositions,
    getLabels,
    getScaleFactor,
    POINT_VISUALIZER_ID,
    updateScatterChart
} from './ThreeUtil';

function getImageUrl(cachedData,
                     categoricalNames,
                     chartOptions,
                     chartSize,
                     markerOpacity,
                     obsCat,
                     pointSize,
                     selection,
                     traceInfo,
                     unselectedMarkerOpacity) {
    let canvas = document.createElement('canvas');
    canvas.width = chartSize * window.devicePixelRatio;
    canvas.height = chartSize * window.devicePixelRatio;
    const context = canvas.getContext('2d');
    context.scale(window.devicePixelRatio, window.devicePixelRatio);
    drawEmbeddingImage(context, {
        width: chartSize,
        height: chartSize
    }, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, chartOptions, categoricalNames, obsCat, cachedData, getSpotRadius(traceInfo, pointSize));
    return canvas.toDataURL();
}

export default function GalleryImage(props) {
    const [url, setUrl] = useState(null);
    const [overlayUrl, setOverlayUrl] = useState(null);
    const [loading, setLoading] = useState(false);
    const elementRef = useRef();


    function onSelect(event) {
        event.preventDefault();
        props.onSelect(props.traceInfo);
    }

    const {
        cachedData,
        categoricalNames,
        chartOptions,
        chartSize,
        containerElement,
        primaryChartSize,
        markerOpacity,
        obsCat,
        pointSize,
        scatterPlot,
        selection,
        traceInfo,
        unselectedMarkerOpacity,
        unselectedPointSize
    } = props;

    useEffect(() => {
        if (traceInfo.type === 'scatter' && traceInfo.embedding.mode == null) {
            let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
            spriteVisualizer.zoomFactor = getScaleFactor(primaryChartSize);

            updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize, unselectedPointSize,
                categoricalNames, chartOptions, obsCat, cachedData, traceInfo.camera);

            const canvas = containerElement.querySelector('canvas');
            const showLabels = obsCat.length > 0 && chartOptions.showGalleryLabels;
            let overlayUrl = null;
            if (showLabels) {
                const labelsPositions = getCategoryLabelsPositions(traceInfo.embedding, obsCat, cachedData);
                const labelCanvas = document.createElement('canvas');
                labelCanvas.width = chartSize * window.devicePixelRatio;
                labelCanvas.height = chartSize * window.devicePixelRatio;
                const context = labelCanvas.getContext('2d');
                context.scale(window.devicePixelRatio, window.devicePixelRatio);
                context.font = 'bold ' + chartOptions.labelFontSize + 'px Roboto Condensed';
                drawLabels(context, getLabels(obsCat, labelsPositions.labels, categoricalNames), labelsPositions.positions, chartOptions, {
                    width: chartSize,
                    height: chartSize
                }, scatterPlot.camera);
                overlayUrl = labelCanvas.toDataURL();
            }
            setUrl(canvas.toDataURL());
            setOverlayUrl(overlayUrl);
            setLoading(false);
        } else if (traceInfo.type === 'image') {
            if (!traceInfo.tileSource.ready) {
                setUrl(null);
                setOverlayUrl(null);
                setLoading(true);
                traceInfo.tileSource.addOnceHandler('ready', () => {
                    setLoading(false);
                    setUrl(getImageUrl(cachedData,
                        categoricalNames,
                        chartOptions,
                        chartSize,
                        markerOpacity,
                        obsCat,
                        pointSize,
                        selection,
                        traceInfo,
                        unselectedMarkerOpacity));
                });
            } else {
                setUrl(getImageUrl(cachedData,
                    categoricalNames,
                    chartOptions,
                    chartSize,
                    markerOpacity,
                    obsCat,
                    pointSize,
                    selection,
                    traceInfo,
                    unselectedMarkerOpacity));
                setOverlayUrl(null);
                setLoading(false);
            }
        } else {
            const containerElement = elementRef.current;
            containerElement.innerHTML = '';
            const svg = traceInfo.gallerySource;
            svg.setAttribute('width', chartSize);
            svg.setAttribute('height', chartSize);
            containerElement.append(svg);
            setUrl(null);
            setOverlayUrl(null);
            setLoading(false);
        }

    }, [containerElement, primaryChartSize, cachedData, categoricalNames, chartOptions, chartSize, markerOpacity, obsCat, pointSize, scatterPlot, selection, traceInfo, unselectedMarkerOpacity, unselectedPointSize]);


    let name = props.traceInfo.name;
    if (name === '__count') {
        name = '';
    }
    return (
        <Box borderColor="text.primary" border={1}
             data-testid="gallery-image"
             style={{display: 'inline-block', margin: 2}}>
            <div style={{
                position: 'relative',
                width: props.chartSize,
                height: props.chartSize,
                cursor: 'pointer'
            }}>
                <Tooltip title={"Embedding: " + props.traceInfo.embedding.name}>
                    <Typography color="textPrimary" variant={"caption"}
                                onClick={onSelect}
                                style={{
                                    marginTop: 3.2,
                                    position: 'absolute',
                                    right: 4,
                                    zIndex: 1000
                                }}>{name}</Typography>
                </Tooltip>
                {loading && <CircularProgress
                    style={{position: 'absolute', left: props.chartSize / 2, top: props.chartSize / 2}}
                    size={20}/>}
                <div onClick={onSelect} ref={elementRef}
                     style={{position: 'absolute', left: 0, top: 0}}></div>
                {url &&
                <div style={{position: 'absolute', left: 0, top: 0}}>
                    <img alt="" src={url}
                         width={props.chartSize * window.devicePixelRatio}
                         height={props.chartSize * window.devicePixelRatio}
                         onClick={onSelect}
                         style={{
                             width: props.chartSize,
                             height: props.chartSize
                         }}/>
                </div>}
                {overlayUrl &&
                <div style={{position: 'absolute', left: 0, top: 0}}>
                    <img alt="" src={overlayUrl}
                         onClick={onSelect}
                         style={{
                             width: props.chartSize,
                             height: props.chartSize
                         }}/></div>}


            </div>
        </Box>
    );


}


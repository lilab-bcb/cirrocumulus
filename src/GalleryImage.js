import {Tooltip} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import CircularProgress from '@material-ui/core/CircularProgress';
import Typography from '@material-ui/core/Typography';
import React from 'react';
import {drawEmbeddingImage, getSpotRadius} from './ImageChart';
import {drawLabels, getVisualizer} from './ScatterChartThree';
import {
    getCategoryLabelsPositions,
    getLabels,
    getScaleFactor,
    POINT_VISUALIZER_ID,
    updateScatterChart
} from './ThreeUtil';


class GalleryImage extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {url: null, overlayUrl: null, loading: false, forceUpdate: false};
        this.elementRef = React.createRef();
    }

    draw() {
        const {
            cachedData,
            categoricalNames,
            chartOptions,
            chartSize,
            containerElement,
            markerOpacity,
            obsCat,
            pointSize,
            scatterPlot,
            selection,
            traceInfo,
            unselectedMarkerOpacity,
        } = this.props;
        if (traceInfo.type === 'scatter') {
            let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
            spriteVisualizer.zoomFactor = this.zoomFactor;
            updateScatterChart(scatterPlot, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, pointSize,
                categoricalNames, chartOptions, obsCat, cachedData);
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

            this.setState({url: canvas.toDataURL(), overlayUrl: overlayUrl, loading: false, element: null});
        } else if (traceInfo.type === 'image') {
            if (!traceInfo.tileSource.ready) {
                this.setState({url: null, overlayUrl: null, loading: true, element: null});
                traceInfo.tileSource.addOnceHandler('ready', () => {
                    this.setState({loading: false});
                });
            } else {
                let canvas = document.createElement('canvas');
                canvas.width = chartSize * window.devicePixelRatio;
                canvas.height = chartSize * window.devicePixelRatio;
                const context = canvas.getContext('2d');
                context.scale(window.devicePixelRatio, window.devicePixelRatio);
                drawEmbeddingImage(context, {
                    width: chartSize,
                    height: chartSize
                }, traceInfo, selection, markerOpacity, unselectedMarkerOpacity, chartOptions, categoricalNames, obsCat, cachedData, getSpotRadius(traceInfo, pointSize));
                this.setState({url: canvas.toDataURL(), overlayUrl: null, loading: false, element: null});
                canvas = null;
            }
        } else {
            const containerElement = this.elementRef.current;
            containerElement.innerHTML = '';
            const svg = traceInfo.gallerySource;
            svg.setAttribute('width', chartSize);
            svg.setAttribute('height', chartSize);
            containerElement.append(svg);
            this.setState({url: null, overlayUrl: null, loading: false});
        }

        // canvas.toBlob(function (blob) {
        //     // let newImg = document.createElement('img');
        //     let url = URL.createObjectURL(blob);
        //     _this.setState({url: url});
        //     // newImg.onload = function () {
        //     //     // no longer need to read the blob so it's revoked
        //     //     URL.revokeObjectURL(url);
        //     // };
        //     //
        //     // newImg.src = url;
        //     // document.body.appendChild(newImg);
        // });

    }


    componentDidMount() {
        this.zoomFactor = getScaleFactor(this.props.primaryChartSize);
        this.draw();
    }


    componentDidUpdate(prevProps, prevState) {
        if (prevProps.primaryChartSize !== this.props.primaryChartSize) {
            this.zoomFactor = getScaleFactor(this.props.primaryChartSize);
        }
        this.draw();
    }


    onSelect = (event) => {
        event.preventDefault();
        this.props.onSelect(this.props.traceInfo);
    };

    render() {

        let name = this.props.traceInfo.name;
        if (name === '__count') {
            name = '';
        }
        return (
            <Box borderColor="text.primary" border={1}
                 style={{display: 'inline-block', margin: 2}}>
                <div style={{
                    position: 'relative',
                    width: this.props.chartSize,
                    height: this.props.chartSize,
                    cursor: 'pointer'
                }}>
                    <Tooltip title={"Embedding: " + this.props.traceInfo.embedding.name}>
                        <Typography color="textPrimary" variant={"caption"}
                                    onClick={this.onSelect}
                                    style={{
                                        marginTop: 3.2,
                                        position: 'absolute',
                                        right: 4,
                                        zIndex: 1000
                                    }}>{name}</Typography>
                    </Tooltip>
                    {this.state.loading && <CircularProgress
                        style={{position: 'absolute', left: this.props.chartSize / 2, top: this.props.chartSize / 2}}
                        size={20}/>}
                    <div onClick={this.onSelect} ref={this.elementRef}
                         style={{position: 'absolute', left: 0, top: 0}}></div>
                    {this.state.url &&
                    <div style={{position: 'absolute', left: 0, top: 0}}>
                        <img alt="" src={this.state.url}
                             width={this.props.chartSize * window.devicePixelRatio}
                             height={this.props.chartSize * window.devicePixelRatio}
                             onClick={this.onSelect}
                             style={{
                                 width: this.props.chartSize,
                                 height: this.props.chartSize
                             }}/>
                    </div>}
                    {this.state.overlayUrl &&
                    <div style={{position: 'absolute', left: 0, top: 0}}>
                        <img alt="" src={this.state.overlayUrl}
                             onClick={this.onSelect}
                             style={{
                                 width: this.props.chartSize,
                                 height: this.props.chartSize
                             }}/></div>}


                </div>
            </Box>
        );

    }
}

export default GalleryImage;


import {Tooltip} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import CircularProgress from '@material-ui/core/CircularProgress';
import Typography from '@material-ui/core/Typography';
import React from 'react';
import {getEmbeddingKey} from './actions';
import {drawImage, getSpotRadius} from './ImageChart';
import {drawLabels, getVisualizer} from './ScatterChartThree';
import {getCategoryLabelsPositions, getScaleFactor, POINT_VISUALIZER_ID, updateScatterChart} from './ThreeUtil';


class GalleryImage extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {url: null, overlayUrl: null, loading: false};
        this.canvasRef = React.createRef();
    }


    draw() {
        const {scatterPlot, chartSize, categoricalNames, containerElement, traceInfo, markerOpacity, unselectedMarkerOpacity, selection, chartOptions, pointSize} = this.props;
        const embedding = traceInfo.embedding;
        const fullName = getEmbeddingKey(embedding);
        const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
        const userPoints = chartSelection ? chartSelection.userPoints : new Set();
        if (!traceInfo.isImage) {
            let spriteVisualizer = getVisualizer(scatterPlot, POINT_VISUALIZER_ID);
            spriteVisualizer.zoomFactor = this.zoomFactor;
            updateScatterChart(scatterPlot, traceInfo, userPoints, markerOpacity, unselectedMarkerOpacity, pointSize,
                categoricalNames, chartOptions);
            const canvas = containerElement.querySelector('canvas');
            const showLabels = chartOptions.showLabels && traceInfo.isCategorical;
            let overlayUrl = null;

            if (showLabels) {
                const labelsPositions = getCategoryLabelsPositions(traceInfo, categoricalNames);
                const labelCanvas = document.createElement('canvas');
                labelCanvas.width = chartSize * window.devicePixelRatio;
                labelCanvas.height = chartSize * window.devicePixelRatio;
                const context = labelCanvas.getContext('2d');
                context.scale(window.devicePixelRatio, window.devicePixelRatio);
                context.font = 'bold ' + chartOptions.labelFontSize + 'px Roboto Condensed';
                drawLabels(context, labelsPositions, chartOptions, {
                    width: chartSize,
                    height: chartSize
                }, scatterPlot.camera);
                overlayUrl = labelCanvas.toDataURL();
            }

            this.setState({url: canvas.toDataURL(), overlayUrl: overlayUrl, loading: false});
        } else {
            if (!traceInfo.tileSource.ready) {
                this.setState({url: null, overlayUrl: null, loading: true});
                traceInfo.tileSource.addOnceHandler('ready', () => {
                    this.setState({loading: false});
                });
            } else {
                let canvas = document.createElement('canvas');
                canvas.width = chartSize * window.devicePixelRatio;
                canvas.height = chartSize * window.devicePixelRatio;
                const context = canvas.getContext('2d');
                context.scale(window.devicePixelRatio, window.devicePixelRatio);
                drawImage(context, {
                    width: chartSize,
                    height: chartSize
                }, traceInfo, userPoints, markerOpacity, unselectedMarkerOpacity, chartOptions, categoricalNames, getSpotRadius(traceInfo, pointSize));
                this.setState({url: canvas.toDataURL(), overlayUrl: null, loading: false});
                canvas = null;
            }
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
            name = 'count';
        }
        return (
            <Box borderColor="text.primary" border={1}
                 style={{display: 'inline-block', margin: 2}}>
                <div style={{position: 'relative', width: this.props.chartSize, height: this.props.chartSize, cursor:'pointer'}}>

                    <Tooltip title={"Embedding: " + this.props.traceInfo.embedding.name}>
                        <Typography color="textPrimary" component={"h4"}
                                    onClick={this.onSelect}
                                    style={{
                                        marginTop: '3.2px',
                                        position: 'absolute',
                                        right: 4,
                                        zIndex: 1000
                                    }}>{name}</Typography>
                    </Tooltip>

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

                    {this.state.loading && <CircularProgress
                        style={{position: 'absolute', left: this.props.chartSize / 2, top: this.props.chartSize / 2}}
                        size={20}/>}
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


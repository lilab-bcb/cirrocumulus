import {Tooltip} from '@material-ui/core';
import Button from '@material-ui/core/Button';
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';
import React from 'react';
import {getEmbeddingKey} from './actions';
import {drawImage} from './ImageChart';
import {updateScatterChart} from './ThreeUtil';


class GalleryImage extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {url: null, forceUpdate: false};
        this.canvasRef = React.createRef();
    }


    draw() {
        let start = new Date().getTime();
        const {scatterPlot, containerElement, traceInfo, markerOpacity, unselectedMarkerOpacity, selection, chartOptions, pointSize} = this.props;
        const embedding = traceInfo.embedding;
        const fullName = getEmbeddingKey(embedding);
        const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
        const userPoints = chartSelection ? chartSelection.userPoints : new Set();
        if (!traceInfo.isImage) {
            updateScatterChart(scatterPlot, traceInfo, userPoints, markerOpacity, unselectedMarkerOpacity, pointSize,
                false, {}, chartOptions.showFog, chartOptions.showAxis, chartOptions.darkMode);
            const canvas = containerElement.querySelector('canvas');
            this.setState({url: canvas.toDataURL()});
        } else {
            if (!traceInfo.tileSource.ready) {
                this.setState({url: null});
                traceInfo.tileSource.addOnceHandler('ready', () => {
                    this.setState((state, props) => {
                        return {forceUpdate: !state.forceUpdate};
                    });
                });
            } else {
                let canvas = document.createElement('canvas');
                canvas.width = this.props.chartSize * window.devicePixelRatio;
                canvas.height = this.props.chartSize * window.devicePixelRatio;
                canvas.style.width = this.props.chartSize * window.devicePixelRatio + ' px';
                canvas.style.height = this.props.chartSize * window.devicePixelRatio + ' px';

                drawImage(canvas.getContext('2d'), {
                    width: this.props.chartSize,
                    height: this.props.chartSize
                }, traceInfo, userPoints, markerOpacity, unselectedMarkerOpacity);
                this.setState({url: canvas.toDataURL()});
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
        this.draw();
    }


    componentDidUpdate(prevProps, prevState) {
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
        return (<Card variant="outlined" style={{display: 'inline-block'}}>
            <CardContent>

                {this.state.url && <img alt="" src={this.state.url}
                                        width={this.props.chartSize * window.devicePixelRatio}
                                        height={this.props.chartSize * window.devicePixelRatio}
                                        style={{
                                            width: this.props.chartSize,
                                            height: this.props.chartSize
                                        }}/>}

            </CardContent>
            <CardActions>
                <Tooltip title={"Embedding: " + this.props.traceInfo.embedding.name}>
                    <Button style={{margin: 'auto'}} onClick={this.onSelect} size="small">{name}</Button>
                </Tooltip>
            </CardActions>
        </Card>);

    }
}

export default GalleryImage;


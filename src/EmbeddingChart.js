import IconButton from '@material-ui/core/IconButton';
import PauseCircleFilledIcon from '@material-ui/icons/PauseCircleFilled';
import PlayCircleFilledIcon from '@material-ui/icons/PlayCircleFilled';
import React from 'react';

import {connect} from 'react-redux';
import {
    getEmbeddingKey,
    handleBrushFilterUpdated,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import createPlotlyComponent from './factory';
import ImageChart from './ImageChart';

const Plot = createPlotlyComponent(window.Plotly);

function xyz2rtz(xyz) {
    return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
    };
}

function rtz2xyz(rtz) {
    return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
    };
}

class EmbeddingChart extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {animating: false, forceRender: true};
    }

    componentWillUnmount() {
        if (this.graphDiv != null) {
            window.Plotly.purge(this.graphDiv);
        }
    }

    onWebglcontextlost() {
        if (this.graphDiv != null) {
            window.Plotly.purge(this.graphDiv);
            window.Plotly.react(this.graphDiv, {
                data: this.props.traceInfo.data,
                layout: this.props.traceInfo.layout,
                config: this.props.config
            });
        }
    }

    toggleAnimation = (event) => {
        this.setState((prevState, props) => ({
            animating: !prevState.animating
        }), () => {
            window.requestAnimationFrame(() => {
                if (this.state.animating) {
                    this.animate();
                }
            });
        });
    };
    animate = () => {
        const eye = this.props.traceInfo.layout.scene.camera.eye;
        const rtz = xyz2rtz(eye);
        rtz.t += Math.PI / 180;
        this.props.traceInfo.layout.scene.camera.eye = rtz2xyz(rtz);
        this.props.traceInfo.layout = Object.assign({}, this.props.traceInfo.layout);
        // force render
        // this.setState((prevState, props) => ({
        //     forceRender: !prevState.forceRender
        // }));
        // hack to update plot because Plot component uses newPlot as re-rendering bug workaround
        window.Plotly.react(this.graphDiv, {
            data: this.props.traceInfo.data,
            layout: this.props.traceInfo.layout,
            config: this.props.config
        });
        if (this.state.animating) {
            window.requestAnimationFrame(() => {
                this.animate();
            });
        }
    };

    onInitialized = (figure, graphDiv) => {
        this.graphDiv = graphDiv;
    };


    onSelect = (event) => {
        const trace = this.props.traceInfo.data[0];
        const name = getEmbeddingKey(trace.embedding);
        if (event == null || event.points == null || event.points.length === 0) {
            this.props.onDeselect({name: name});
            return;
        }
        let userPath;
        if (event.range) { // rect
            userPath = {
                shape: 'rect',
                x: event.range.x[0],
                y: event.range.y[0],
                width: event.range.x[1] - event.range.x[0],
                height: event.range.y[1] - event.range.y[0]
            };
        } else { // lasso
            userPath = [];
            for (let i = 0; i < event.lassoPoints.x.length; i++) {
                userPath.push([event.lassoPoints.x[i], event.lassoPoints.y[i]]);
            }
        }
        this.props.onSelect({
            name: name,
            value: {basis: trace.embedding, selectedpoints: event.points, path: userPath}
        });
    };

    onDeselect = (event) => {
        const trace = this.props.traceInfo.data[0];
        const name = getEmbeddingKey(trace.embedding);
        this.props.onDeselect({name: name});

    };

    render() {
        const {traceInfo, config, shape, nObsSelected, onDeselect, onSelect, globalFeatureSummary, featureSummary, datasetFilter, handleColorChange, handleDimensionFilterUpdated, handleMeasureFilterUpdated} = this.props;

        return (
            <div style={this.props.style}>
                {!traceInfo.data[0].isImage && <Plot
                    style={{display: 'inline-block'}}
                    data={traceInfo.data}
                    onInitialized={this.onInitialized}
                    layout={traceInfo.layout}
                    config={config}
                    onDeselect={this.onDeselect}
                    onWebglcontextlost={this.onWebglcontextlost}
                    onSelected={this.onSelect}
                />}


                {/*{!traceInfo.data[0].isImage && traceInfo.data[0].type === 'scattergl' && <Scatter2d*/}
                {/*    style={{display: 'inline-block'}}*/}
                {/*    data={traceInfo.data}*/}
                {/*    onInitialized={this.onInitialized}*/}
                {/*    layout={traceInfo.layout}*/}
                {/*    config={config}*/}
                {/*    onDeselect={onDeselect}*/}
                {/*    onSelected={onSelect}*/}
                {/*/>}*/}

                {traceInfo.data[0].isImage && <ImageChart
                    style={{display: 'inline-block'}}
                    data={traceInfo.data}
                    onInitialized={this.onInitialized}
                    layout={traceInfo.layout}
                    config={config}
                    onDeselect={onDeselect}
                    onSelected={onSelect}
                />}
                {traceInfo.data[0].type === 'scatter3d' && !this.state.animating &&
                <IconButton onClick={this.toggleAnimation} aria-label="Play">
                    <PlayCircleFilledIcon/>
                </IconButton>}

                {traceInfo.data[0].type === 'scatter3d' && this.state.animating &&
                <IconButton onClick={this.toggleAnimation} aria-label="Pause">
                    <PauseCircleFilledIcon/>
                </IconButton>}


                {traceInfo.continuous ?
                    <ColorSchemeLegendWrapper
                        width={140}
                        label={true}
                        height={30}
                        handleUpdate={handleMeasureFilterUpdated}
                        datasetFilter={datasetFilter}
                        scale={traceInfo.colorScale}
                        featureSummary={featureSummary}
                        globalFeatureSummary={globalFeatureSummary}
                        nObs={shape[0]}
                        nObsSelected={nObsSelected}
                        maxHeight={traceInfo.layout.height}
                        name={traceInfo.name}
                    /> :
                    <CategoricalLegend datasetFilter={datasetFilter}
                                       handleClick={handleDimensionFilterUpdated}
                                       handleColorChange={handleColorChange}
                                       name={traceInfo.name}
                                       scale={traceInfo.colorScale}
                                       maxHeight={traceInfo.layout.height - 24}
                                       clickEnabled={true}
                                       nObs={shape[0]}
                                       nObsSelected={nObsSelected}
                                       globalFeatureSummary={globalFeatureSummary}
                                       featureSummary={featureSummary}/>}

            </div>);

    }
}

const mapStateToProps = state => {
    return {
        numberOfBins: state.numberOfBins,
        binValues: state.binValues,
        embeddingChartSize: state.embeddingChartSize,
        config: state.plotConfig,
        datasetFilter: state.datasetFilter,
        featureSummary: state.featureSummary,
        shape: state.dataset.shape,
        nObsSelected: state.selection.count,
        globalFeatureSummary: state.globalFeatureSummary
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handleDimensionFilterUpdated: (e) => {
            dispatch(handleDimensionFilterUpdated(e));
        },
        handleColorChange: (e) => {
            dispatch(handleColorChange(e));
        },
        handleMeasureFilterUpdated: (e) => {
            dispatch(handleMeasureFilterUpdated(e));
        },
        onSelect: (e) => {
            console.log(e);

            // name is full basis name
            dispatch(handleBrushFilterUpdated(e));
        },
        onDeselect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingChart));


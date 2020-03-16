import React from 'react';


import {connect} from 'react-redux';
import {
    handleBrushFilterUpdated,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import {getChartSize} from './PlotUtil';
import ScatterChartThree from './ScatterChartThree';


class EmbeddingChart extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {animating: false};
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


    render() {
        const {
            traceInfo, config, shape, nObsSelected, onDeselect, onSelect, globalFeatureSummary, featureSummary,
            datasetFilter, handleColorChange, handleDimensionFilterUpdated, handleMeasureFilterUpdated
        } = this.props;

        const chartSize = getChartSize();
        return (

            <div style={this.props.style}>
                <div style={{position: 'absolute', right: 10}}>
                    {traceInfo.continuous ?
                        <ColorSchemeLegendWrapper
                            width={140}
                            label={true}
                            showColorScheme={true}
                            height={30}
                            handleUpdate={handleMeasureFilterUpdated}
                            datasetFilter={datasetFilter}
                            scale={traceInfo.colorScale}
                            featureSummary={featureSummary}
                            globalFeatureSummary={globalFeatureSummary}
                            nObs={shape[0]}
                            nObsSelected={nObsSelected}
                            name={traceInfo.name}
                        /> :
                        <CategoricalLegend datasetFilter={datasetFilter}
                                           handleClick={handleDimensionFilterUpdated}
                                           handleColorChange={handleColorChange}
                                           name={traceInfo.name}
                                           scale={traceInfo.colorScale}
                                           maxHeight={chartSize.height - 24}
                                           clickEnabled={true}
                                           nObs={shape[0]}
                                           nObsSelected={nObsSelected}
                                           globalFeatureSummary={globalFeatureSummary}
                                           featureSummary={featureSummary}/>}
                </div>

                {/*{!traceInfo.data[0].isImage && <Plot*/}
                {/*    style={{display: 'inline-block'}}*/}
                {/*    data={traceInfo.data}*/}
                {/*    onInitialized={this.onInitialized}*/}
                {/*    layout={traceInfo.layout}*/}
                {/*    config={config}*/}
                {/*    onDeselect={this.onDeselect}*/}
                {/*    onWebglcontextlost={this.onWebglcontextlost}*/}
                {/*    onSelected={this.onSelect}*/}
                {/*/>}*/}

                {!traceInfo.isImage &&
                <ScatterChartThree traceInfo={traceInfo}
                                   selection={this.props.selection}
                                   onDeselect={this.props.onDeselect}
                                   onSelected={this.props.onSelect}
                                   markerOpacity={this.props.markerOpacity}
                                   color={traceInfo.marker.color}
                                   unselectedMarkerOpacity={this.props.unselectedMarkerOpacity}
                />}


                {traceInfo.isImage && <ImageChart
                    style={{display: 'inline-block'}}
                    data={traceInfo.data}
                    onInitialized={this.onInitialized}
                    layout={traceInfo.layout}
                    config={config}
                    onDeselect={onDeselect}
                    onSelected={onSelect}
                />}
                {/*{traceInfo.data[0].type === 'scatter3d' && !this.state.animating &&*/}
                {/*<IconButton onClick={this.toggleAnimation} aria-label="Play">*/}
                {/*    <PlayCircleFilledIcon/>*/}
                {/*</IconButton>}*/}

                {/*{traceInfo.data[0].type === 'scatter3d' && this.state.animating &&*/}
                {/*<IconButton onClick={this.toggleAnimation} aria-label="Pause">*/}
                {/*    <PauseCircleFilledIcon/>*/}
                {/*</IconButton>}*/}


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


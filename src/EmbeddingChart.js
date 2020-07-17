import React from 'react';


import {connect} from 'react-redux';
import {
    handleBrushFilterUpdated,
    handleCategoricalNameChange,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import ScatterChartThree from './ScatterChartThree';
import {getChartSize} from './util';


class EmbeddingChart extends React.PureComponent {


    render() {
        const {
            traceInfo, categoricalNames, config, shape, nObsSelected, onDeselect, onSelect, globalFeatureSummary, featureSummary,
            datasetFilter, handleColorChange, handleNameChange, handleDimensionFilterUpdated, handleMeasureFilterUpdated
        } = this.props;

        const chartSize = getChartSize();
        return (

            <div style={{position: 'relative'}}>
                <div style={{position: 'absolute', right: 10, zIndex: 1000}}>
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
                                           handleNameChange={handleNameChange}
                                           name={traceInfo.name}
                                           categoricalNames={categoricalNames}
                                           scale={traceInfo.colorScale}
                                           maxHeight={chartSize.height - 24}
                                           clickEnabled={true}
                                           nObs={shape[0]}
                                           nObsSelected={nObsSelected}
                                           globalFeatureSummary={globalFeatureSummary}
                                           featureSummary={featureSummary}/>}
                </div>


                {!traceInfo.isImage &&
                <ScatterChartThree traceInfo={traceInfo}
                                   categoricalNames={categoricalNames}
                                   selection={this.props.selection}
                                   onDeselect={this.props.onDeselect}
                                   onSelected={this.props.onSelect}
                                   pointSize={this.props.pointSize}
                                   markerOpacity={this.props.markerOpacity}
                                   color={traceInfo.colors}
                                   onGallery={this.props.onGallery}
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


            </div>);

    }
}

const mapStateToProps = state => {
    return {
        categoricalNames: state.categoricalNames,
        numberOfBins: state.numberOfBins,
        binValues: state.binValues,
        embeddingChartSize: state.embeddingChartSize,
        config: state.plotConfig,
        datasetFilter: state.datasetFilter,
        featureSummary: state.featureSummary,
        shape: state.dataset.shape,
        nObsSelected: state.selection.count,
        pointSize: state.pointSize,
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
        handleNameChange: (e) => {
            dispatch(handleCategoricalNameChange(e));
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


import Paper from '@material-ui/core/Paper';
import React from 'react';


import {connect} from 'react-redux';
import {
    handleBrushFilterUpdated,
    handleCategoricalNameChange,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleMeasureFilterUpdated,
    setChartOptions
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import ScatterChartThree from './ScatterChartThree';
import {getChartSize} from './util';


class EmbeddingChart extends React.PureComponent {


    render() {
        const {
            traceInfo, categoricalNames, shape, nObsSelected, globalFeatureSummary, featureSummary,
            datasetFilter, handleColorChange, handleNameChange, handleDimensionFilterUpdated, handleMeasureFilterUpdated
        } = this.props;
        const chartSize = getChartSize();
        return (

            <div style={{position: 'relative'}}>
                <Paper elevation={0} style={{position: 'absolute', right: 10, zIndex: 1000}}>
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
                </Paper>


                {!traceInfo.isImage &&
                <ScatterChartThree traceInfo={traceInfo}
                                   setChartOptions={this.props.handleChartOptions}
                                   chartOptions={this.props.chartOptions}
                                   categoricalNames={categoricalNames}
                                   selection={this.props.selection}
                                   onDeselect={this.props.onDeselect}
                                   onSelected={this.props.onSelect}
                                   pointSize={this.props.pointSize}
                                   markerOpacity={this.props.markerOpacity}
                                   unselectedMarkerOpacity={this.props.unselectedMarkerOpacity}
                                   color={traceInfo.colors}
                                   onGallery={this.props.onGallery}

                />}


                {traceInfo.isImage && <ImageChart
                    setChartOptions={this.props.handleChartOptions}
                    chartOptions={this.props.chartOptions}
                    style={{display: 'inline-block'}}
                    traceInfo={traceInfo}
                    categoricalNames={categoricalNames}
                    selection={this.props.selection}
                    onInitialized={this.onInitialized}
                    markerOpacity={this.props.markerOpacity}
                    unselectedMarkerOpacity={this.props.unselectedMarkerOpacity}
                    onGallery={this.props.onGallery}
                    onDeselect={this.props.onDeselect}
                    onSelected={this.props.onSelect}
                />}


            </div>);

    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        categoricalNames: state.categoricalNames,
        numberOfBins: state.numberOfBins,
        binValues: state.binValues,
        embeddingChartSize: state.embeddingChartSize,
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
        handleChartOptions: (options) => {
            dispatch(setChartOptions(options));
        },
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


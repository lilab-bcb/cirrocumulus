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
import ImageChart from './ImageChart';
import ScatterChartThree from './ScatterChartThree';


class EmbeddingChart extends React.PureComponent {

    setChartRef = (element) => {
        this.props.chartOptions.ref = element;
    };


    render() {
        const {traceInfo, categoricalNames, primaryChartSize} = this.props;
        const traceName = traceInfo.name === '__count' ? '' : traceInfo.name;
        return (

            <div style={{position: 'relative'}}>
                <Paper elevation={0} style={{position: 'absolute', right: 10, zIndex: 1000}}>
                    <h4 style={{marginTop: '3.2px'}}>{traceName}</h4>
                    {/*{traceInfo.continuous ?*/}
                    {/*    <ColorSchemeLegendWrapper*/}
                    {/*        width={140}*/}
                    {/*        label={true}*/}
                    {/*        showColorScheme={true}*/}
                    {/*        height={30}*/}
                    {/*        handleDomain={handleDomain}*/}
                    {/*        handleUpdate={handleMeasureFilterUpdated}*/}
                    {/*        datasetFilter={datasetFilter}*/}
                    {/*        scale={traceInfo.colorScale}*/}
                    {/*        featureSummary={featureSummary}*/}
                    {/*        globalFeatureSummary={globalFeatureSummary}*/}
                    {/*        nObs={shape[0]}*/}
                    {/*        nObsSelected={nObsSelected}*/}
                    {/*        name={traceInfo.name}*/}
                    {/*    /> :*/}
                    {/*    <CategoricalLegend datasetFilter={datasetFilter}*/}
                    {/*                       handleClick={handleDimensionFilterUpdated}*/}
                    {/*                       handleColorChange={handleColorChange}*/}
                    {/*                       handleNameChange={handleNameChange}*/}
                    {/*                       name={traceInfo.name}*/}
                    {/*                       categoricalNames={categoricalNames}*/}
                    {/*                       scale={traceInfo.colorScale}*/}
                    {/*                       maxHeight={primaryChartSize.height - 24}*/}
                    {/*                       clickEnabled={true}*/}
                    {/*                       nObs={shape[0]}*/}
                    {/*                       nObsSelected={nObsSelected}*/}
                    {/*                       globalFeatureSummary={globalFeatureSummary}*/}
                    {/*                       featureSummary={featureSummary}/>}*/}
                </Paper>


                {!traceInfo.isImage &&
                <ScatterChartThree traceInfo={traceInfo}
                                   chartSize={primaryChartSize}
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
                                   onMoreOptions={this.props.onMoreOptions}
                                   ref={this.setChartRef}

                />}


                {traceInfo.isImage && <ImageChart
                    setChartOptions={this.props.handleChartOptions}
                    chartOptions={this.props.chartOptions}
                    style={{display: 'inline-block'}}
                    traceInfo={traceInfo}
                    pointSize={this.props.pointSize}
                    chartSize={primaryChartSize}
                    categoricalNames={categoricalNames}
                    selection={this.props.selection}
                    onInitialized={this.onInitialized}
                    markerOpacity={this.props.markerOpacity}
                    unselectedMarkerOpacity={this.props.unselectedMarkerOpacity}
                    onDeselect={this.props.onDeselect}
                    onSelected={this.props.onSelect}
                    onGallery={this.props.onGallery}
                    onMoreOptions={this.props.onMoreOptions}
                    ref={this.setChartRef}
                />}


            </div>);

    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        primaryChartSize: state.primaryChartSize,
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


import {Tooltip, Typography} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Link from '@material-ui/core/Link';
import React from 'react';


import {connect} from 'react-redux';
import {
    handleBrushFilterUpdated,
    handleCategoricalNameChange,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleDomainChange,
    handleMeasureFilterUpdated,
    MORE_OPTIONS_DIALOG,
    setChartOptions,
    setDialog,
    setPrimaryChartSize
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import MetaEmbedding from './MetaEmbedding';
import ScatterChartThree from './ScatterChartThree';
import {splitSearchTokens} from './util';

class EmbeddingChart extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {showDetails: true};
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.traceInfo.name !== this.props.traceInfo.name) {
            this.setState({showDetails: true});
        }
    }

    handleExpandClick = (e) => {
        e.preventDefault();
        this.setState({showDetails: !this.state.showDetails});
    };

    render() {
        const {
            cachedData,
            categoricalNames,
            chartOptions,
            dataset,
            datasetFilter,
            embeddingLabels,
            featureSummary,
            globalFeatureSummary,
            markerOpacity,
            nObsSelected,
            onChartOptions,
            onColorChange,
            onDeselect,
            onDimensionFilterUpdated,
            onDomain,
            onGallery,
            onMeasureFilterUpdated,
            onNameChange,
            onSelect,
            pointSize,
            primaryChartSize,
            searchTokens,
            selection,
            shape,
            traceInfo,
            unselectedMarkerOpacity
        } = this.props;

        const activeEmbeddingLabels = splitSearchTokens(searchTokens).obsCat.filter(item => embeddingLabels.indexOf(item) !== -1);
        const displayName = traceInfo.name === '__count' ? '' : traceInfo.name;

        return (
            <div style={{position: 'relative'}}>
                <Box color="text.primary" style={{
                    marginTop: '3.2px',
                    position: 'absolute',
                    textAlign: 'right',
                    right: 8,
                    zIndex: 1000
                }}>
                    {displayName !== '' ? <Tooltip title={"Embedding: " + traceInfo.embedding.name}>
                        <Link onClick={this.handleExpandClick}>
                            <Typography
                                color="textPrimary" style={{marginRight: 14}}
                                component={"h4"}>{displayName} {!traceInfo.continuous ?
                                <small>({globalFeatureSummary[traceInfo.name].categories.length})</small> : null}</Typography></Link>
                    </Tooltip> : null}

                    {traceInfo.continuous ?
                        <ColorSchemeLegendWrapper
                            key={traceInfo.name}
                            handleDomain={onDomain}
                            width={140}
                            showColorScheme={false}
                            height={30}
                            style={{
                                display: this.state.showDetails ? 'block' : 'none',
                            }}
                            handleUpdate={onMeasureFilterUpdated}
                            datasetFilter={datasetFilter}
                            scale={traceInfo.colorScale}
                            featureSummary={featureSummary}
                            globalFeatureSummary={globalFeatureSummary}
                            nObs={shape[0]}
                            nObsSelected={nObsSelected}
                            maxHeight={null}
                            name={traceInfo.name}
                        /> :
                        <CategoricalLegend
                            key={traceInfo.name}
                            style={{
                                display: this.state.showDetails ? 'block' : 'none',
                            }}
                            datasetFilter={datasetFilter}
                            handleClick={onDimensionFilterUpdated}
                            handleColorChange={onColorChange}
                            handleNameChange={onNameChange}
                            categoricalNames={categoricalNames}
                            name={traceInfo.name}
                            scale={traceInfo.colorScale}
                            maxHeight={300}
                            clickEnabled={true}
                            nObs={shape[0]}
                            nObsSelected={nObsSelected}
                            globalFeatureSummary={globalFeatureSummary}
                            featureSummary={featureSummary}/>
                    }
                </Box>


                {traceInfo.type === 'scatter' &&
                <ScatterChartThree traceInfo={traceInfo}
                                   cachedData={cachedData}
                                   obsCat={activeEmbeddingLabels}
                                   chartSize={primaryChartSize}
                                   setChartOptions={onChartOptions}
                                   chartOptions={chartOptions}
                                   categoricalNames={categoricalNames}
                                   selection={selection}
                                   onDeselect={onDeselect}
                                   onSelected={onSelect}
                                   pointSize={pointSize}
                                   markerOpacity={markerOpacity}
                                   unselectedMarkerOpacity={unselectedMarkerOpacity}
                                   color={traceInfo.colors}
                                   onGallery={onGallery}

                />}

                {traceInfo.type === 'meta_image' &&
                <MetaEmbedding traceInfo={traceInfo}
                               cachedData={cachedData}
                               chartSize={primaryChartSize}
                               setChartOptions={onChartOptions}
                               chartOptions={chartOptions}
                               dataset={dataset}
                               categoricalNames={categoricalNames}
                               markerOpacity={markerOpacity}
                               onGallery={onGallery}

                />}
                {traceInfo.type === 'image' && <ImageChart
                    cachedData={cachedData}
                    obsCat={activeEmbeddingLabels}
                    setChartOptions={onChartOptions}
                    chartOptions={chartOptions}
                    style={{display: 'inline-block'}}
                    traceInfo={traceInfo}
                    pointSize={pointSize}
                    chartSize={primaryChartSize}
                    categoricalNames={categoricalNames}
                    selection={selection}
                    onInitialized={this.onInitialized}
                    markerOpacity={markerOpacity}
                    unselectedMarkerOpacity={unselectedMarkerOpacity}
                    onDeselect={onDeselect}
                    onSelected={onSelect}
                    onGallery={onGallery}
                />}


            </div>);

    }
}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        dataset: state.dataset,
        datasetFilter: state.datasetFilter,
        embeddingChartSize: state.embeddingChartSize,
        embeddingData: state.embeddingData,
        embeddingLabels: state.embeddingLabels,
        featureSummary: state.featureSummary,
        globalFeatureSummary: state.globalFeatureSummary,
        markerOpacity: state.markerOpacity,
        nObsSelected: state.selection.count,
        pointSize: state.pointSize,
        primaryChartSize: state.primaryChartSize,
        primaryTraceKey: state.primaryTraceKey,
        shape: state.dataset.shape,
        searchTokens: state.searchTokens,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onChartOptions: (options) => {
            dispatch(setChartOptions(options));
        },
        onDomain: (value) => {
            dispatch(handleDomainChange(value));
        },
        onDimensionFilterUpdated: (e) => {
            dispatch(handleDimensionFilterUpdated(e));
        },
        onColorChange: (e) => {
            dispatch(handleColorChange(e));
        },
        onNameChange: (e) => {
            dispatch(handleCategoricalNameChange(e));
        },
        onMeasureFilterUpdated: (e) => {
            dispatch(handleMeasureFilterUpdated(e));
        },
        onSelect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        },
        onDeselect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        },

    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingChart));


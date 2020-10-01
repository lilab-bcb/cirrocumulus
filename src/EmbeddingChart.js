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
import ScatterChartThree from './ScatterChartThree';

class EmbeddingChart extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {showDetails: true};
    }

    setChartRef = (element) => {
        this.props.chartOptions.ref = element;

    };

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
            onChartOptions, onMoreOptions, onDomain, onDimensionFilterUpdated,
            onColorChange, onNameChange, onMeasureFilterUpdated, onSelect, onDeselect, onGallery,
            traceInfo, selection, datasetFilter, chartOptions, featureSummary, markerOpacity, unselectedMarkerOpacity, pointSize, globalFeatureSummary, shape, nObsSelected, categoricalNames, primaryChartSize
        } = this.props;
        const displayName = traceInfo.name === '__count' ? '' : traceInfo.name;

        return (
            <div  style={{position: 'relative'}}>
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


                {!traceInfo.isImage &&
                <ScatterChartThree traceInfo={traceInfo}
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
                                   onMoreOptions={onMoreOptions}
                                   ref={this.setChartRef}

                />}


                {traceInfo.isImage && <ImageChart
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
                    onMoreOptions={onMoreOptions}
                    ref={this.setChartRef}
                />}


            </div>);

    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        embeddingChartSize: state.embeddingChartSize,
        pointSize: state.pointSize,
        embeddingData: state.embeddingData,
        primaryChartSize: state.primaryChartSize,
        markerOpacity: state.markerOpacity,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        primaryTraceKey: state.primaryTraceKey,
        categoricalNames: state.categoricalNames,
        globalFeatureSummary: state.globalFeatureSummary,
        featureSummary: state.featureSummary,
        shape: state.dataset.shape,
        nObsSelected: state.selection.count,
        dataset: state.dataset,
        datasetFilter: state.datasetFilter
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onChartOptions: (options) => {
            dispatch(setChartOptions(options));
        },
        onMoreOptions: () => {
            dispatch(setDialog(MORE_OPTIONS_DIALOG));
        },
        onPrimaryChartSize: value => {
            dispatch(setPrimaryChartSize(value));
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


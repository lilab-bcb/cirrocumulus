import {Tooltip, Typography} from '@material-ui/core';
import Box from '@material-ui/core/Box';
import Link from '@material-ui/core/Link';
import {find} from 'lodash';
import React from 'react';


import {connect} from 'react-redux';
import {
    getTraceKey,
    handleBrushFilterUpdated,
    handleCategoricalNameChange,
    handleColorChange,
    handleDimensionFilterUpdated,
    handleDomainChange,
    handleMeasureFilterUpdated,
    setChartOptions,
    setEmbeddingData,
    setPrimaryChartSize
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import MetaEmbedding from './MetaEmbedding';
import ScatterChartThree from './ScatterChartThree';
import {splitSearchTokens, TRACE_TYPE_META_IMAGE} from './util';
import memoize from 'memoize-one';

// TODO-this causes an unnecessary redraw when obsCat is updated
const getActiveEmbeddingLabels = memoize(
    (searchTokens, embeddingLabels) => {
        return splitSearchTokens(searchTokens).obsCat.filter(item => embeddingLabels.indexOf(item) !== -1);
    }
);

class EmbeddingChart extends React.PureComponent {


    constructor(props) {
        super(props);
        this.state = {showLegend: true};
        this.resizeListener = () => {
            let width = window.innerWidth - 280;
            let height = Math.max(300, window.innerHeight - 370);
            this.props.handlePrimaryChartSize({width: width, height: height});
            this.windowHeight = window.innerHeight;
        };
        window.addEventListener('resize', this.resizeListener);
    }

    componentWillUnmount() {
        window.removeEventListener('resize', this.resizeListener);
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (prevProps.activeFeature == null || this.props.activeFeature == null || prevProps.activeFeature.name !== this.props.activeFeature.name) {
            this.setState({showLegend: true});
        }
    }

    handleToggleLegend = (e) => {
        e.preventDefault();
        this.setState({showLegend: !this.state.showLegend});
    };

    onCamera = (eventName, cameraDef) => {
        const {embeddingData, activeFeature} = this.props;
        const primaryTrace = find(embeddingData, item => getTraceKey(item) === activeFeature.embeddingKey);
        for (let i = 0, n = embeddingData.length; i < n; i++) {
            if (primaryTrace.embedding.name === embeddingData[i].embedding.name) {
                embeddingData[i].camera = cameraDef;
            }
        }
        this.props.handleEmbeddingData(this.props.embeddingData.slice());
    };


    render() {
        const {
            activeFeature,
            cachedData,
            categoricalNames,
            chartOptions,
            dataset,
            datasetFilter,
            embeddingLabels,
            embeddingData,
            featureSummary,
            globalFeatureSummary,
            markerOpacity,
            onChartOptions,
            onColorChange,
            onDimensionFilterUpdated,
            onDomain,
            onGallery,
            onMeasureFilterUpdated,
            onCategoricalNameChange,
            onSelect,
            pointSize,
            primaryChartSize,
            searchTokens,
            selection,
            setTooltip,
            shape,
            unselectedMarkerOpacity,
            unselectedPointSize
        } = this.props;

        if (activeFeature == null) {
            return null;
        }
        const primaryTrace = find(embeddingData, item => getTraceKey(item) === activeFeature.embeddingKey);
        if (primaryTrace == null) {
            return null;
        }
        const nObsSelected = selection.size;
        const activeEmbeddingLabels = getActiveEmbeddingLabels(searchTokens, embeddingLabels);
        const displayName = primaryTrace.name === '__count' ? '' : primaryTrace.name;

        return (
            <div style={{position: 'relative'}}>
                <Box data-testid="chart-extra" color="text.primary" style={{
                    marginTop: 3.2,
                    position: 'absolute',
                    textAlign: 'right',
                    overflow: 'hidden',
                    // whiteSpace: 'nowrap',
                    textOverflow: 'ellipsis',
                    maxWidth: 300,
                    right: 8,
                    zIndex: 1000
                }}>
                    {displayName !== '' &&
                    <Tooltip title={"Embedding: " + primaryTrace.embedding.name}><Link
                        onClick={this.handleToggleLegend}>
                        <Typography
                            color="textPrimary" style={{marginRight: 14}}
                            component={"h4"}>{displayName}</Typography></Link></Tooltip>
                    }
                    {primaryTrace.continuous ?
                        <ColorSchemeLegendWrapper
                            key={primaryTrace.name}
                            handleDomain={onDomain}
                            width={140}
                            showColorScheme={false}
                            height={30}
                            style={{
                                display: this.state.showLegend ? 'block' : 'none'
                            }}
                            handleUpdate={onMeasureFilterUpdated}
                            datasetFilter={datasetFilter}
                            colorScale={primaryTrace.colorScale}
                            featureSummary={featureSummary}
                            globalFeatureSummary={globalFeatureSummary}
                            nObs={shape[0]}
                            nObsSelected={nObsSelected}
                            maxHeight={null}
                            name={primaryTrace.name}
                        /> :
                        <CategoricalLegend
                            key={primaryTrace.name}
                            visible={this.state.showLegend}
                            height={primaryChartSize.height - 40}
                            features={dataset.features}
                            datasetFilter={datasetFilter}
                            handleClick={onDimensionFilterUpdated}
                            handleColorChange={onColorChange}
                            handleNameChange={onCategoricalNameChange}
                            categoricalNames={categoricalNames}
                            name={primaryTrace.name}
                            scale={primaryTrace.colorScale}
                            nObs={shape[0]}
                            nObsSelected={nObsSelected}
                            globalFeatureSummary={globalFeatureSummary}
                            featureSummary={featureSummary}/>
                    }
                </Box>

                {primaryTrace.type === 'scatter' &&
                <ScatterChartThree trace={primaryTrace}
                                   cachedData={cachedData}
                                   obsCat={activeEmbeddingLabels}
                                   chartSize={primaryChartSize}
                                   setChartOptions={onChartOptions}
                                   chartOptions={chartOptions}
                                   categoricalNames={categoricalNames}
                                   selection={selection}
                                   onSelected={onSelect}
                                   pointSize={pointSize}
                                   unselectedPointSize={unselectedPointSize}
                                   markerOpacity={markerOpacity}
                                   unselectedMarkerOpacity={unselectedMarkerOpacity}
                                   color={primaryTrace.colors}
                                   onGallery={onGallery}
                                   onCamera={this.onCamera}
                                   setTooltip={setTooltip}

                />}

                {primaryTrace.type === TRACE_TYPE_META_IMAGE &&
                <MetaEmbedding trace={primaryTrace}
                               cachedData={cachedData}
                               chartSize={primaryChartSize}
                               setChartOptions={onChartOptions}
                               chartOptions={chartOptions}
                               dataset={dataset}
                               selection={selection}
                               categoricalNames={categoricalNames}
                               markerOpacity={markerOpacity}
                               onSelected={onSelect}
                               onGallery={onGallery}
                               setTooltip={setTooltip}

                />}
                {primaryTrace.type === 'image' && <ImageChart
                    cachedData={cachedData}
                    obsCat={activeEmbeddingLabels}
                    setChartOptions={onChartOptions}
                    chartOptions={chartOptions}
                    style={{display: 'inline-block'}}
                    trace={primaryTrace}
                    pointSize={pointSize}
                    chartSize={primaryChartSize}
                    categoricalNames={categoricalNames}
                    selection={selection}
                    onInitialized={this.onInitialized}
                    markerOpacity={markerOpacity}
                    unselectedMarkerOpacity={unselectedMarkerOpacity}
                    onSelected={onSelect}
                    onGallery={onGallery}
                    setTooltip={setTooltip}
                />}
            </div>);

    }
}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        embeddingData: state.embeddingData,
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        dataset: state.dataset,
        datasetFilter: state.datasetFilter,
        embeddingLabels: state.embeddingLabels,
        featureSummary: state.featureSummary,
        globalFeatureSummary: state.globalFeatureSummary,
        markerOpacity: state.markerOpacity,
        selection: state.selection,
        pointSize: state.pointSize,
        primaryChartSize: state.primaryChartSize,
        shape: state.dataset.shape,
        searchTokens: state.searchTokens,
        unselectedMarkerOpacity: state.unselectedMarkerOpacity,
        unselectedPointSize: state.unselectedPointSize
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
        onCategoricalNameChange: (e) => {
            dispatch(handleCategoricalNameChange(e));
        },
        onMeasureFilterUpdated: (e) => {
            dispatch(handleMeasureFilterUpdated(e));
        },
        onSelect: (e) => {
            dispatch(handleBrushFilterUpdated(e));
        },
        handlePrimaryChartSize: value => {
            dispatch(setPrimaryChartSize(value));
        },
        handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(EmbeddingChart));


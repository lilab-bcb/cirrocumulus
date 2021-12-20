import {Tooltip, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import {find} from 'lodash';
import React, {useEffect, useState} from 'react';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';

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
    setLegendScrollPosition,
    setSearchTokens,
    setWindowSize
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import MetaEmbedding from './MetaEmbedding';
import ScatterChartThree from './ScatterChartThree';
import {FEATURE_TYPE, TRACE_TYPE_META_IMAGE} from './util';
import memoize from 'memoize-one';


const getActiveEmbeddingLabels = memoize(
    (searchTokens, embeddingLabels) => {
        return searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT && embeddingLabels.indexOf(item.id) !== -1).map(item => item.id);
    }
);

function EmbeddingChart(props) {
    const {
        activeFeature,
        cachedData,
        categoricalNames,
        chartOptions,
        dataset,
        datasetFilter,
        embeddingData,
        embeddingLabels,
        featureSummary,
        globalFeatureSummary,
        handleEmbeddingData,
        handleScrollPosition,
        handleSearchTokens,
        handleWindowSize,
        legendScrollPosition,
        markerOpacity,
        onCategoricalNameChange,
        onChartOptions,
        onColorChange,
        onDimensionFilterUpdated,
        onDomain,
        onGallery,
        onMeasureFilterUpdated,
        onSelect,
        pointSize,
        primaryChartSize,
        searchTokens,
        selection,
        serverInfo,
        shape,
        unselectedMarkerOpacity,
        unselectedPointSize
    } = props;

    const [showLegend, setShowLegend] = useState(true);
    
    useEffect(() => {
        window.addEventListener('resize', handleWindowSize);
        return () => {
            window.removeEventListener('resize', handleWindowSize);
        };
    }, [handleWindowSize]);


    function handleToggleLegend(e) {
        e.preventDefault();
        setShowLegend(!showLegend);
    }

    function onAddFeatures(features) {
        const values = searchTokens.filter(token => token.type !== FEATURE_TYPE.X);
        const xTokens = searchTokens.filter(token => token.type === FEATURE_TYPE.X);
        const existingXFeatures = new Set();
        xTokens.forEach(item => existingXFeatures.add(item.id));
        features.forEach(feature => {
            if (!existingXFeatures.has(feature)) {
                xTokens.push({id: feature, type: FEATURE_TYPE.X});
                existingXFeatures.add(feature);
            }
        });
        handleSearchTokens(values.concat(xTokens));
    }

    function onCamera(eventName, cameraDef) {
        const primaryTrace = find(embeddingData, item => getTraceKey(item) === activeFeature.embeddingKey);
        for (let i = 0, n = embeddingData.length; i < n; i++) {
            if (primaryTrace.embedding.name === embeddingData[i].embedding.name) {
                embeddingData[i].camera = cameraDef;
            }
        }
        handleEmbeddingData(embeddingData.slice());
    }


    if (activeFeature == null) {
        return null;
    }
    const primaryTrace = find(embeddingData, item => getTraceKey(item) === activeFeature.embeddingKey);
    if (primaryTrace == null) {
        console.log(activeFeature.embeddingKey + ' not found');
        return null;
    }
    const nObsSelected = selection != null ? selection.size : 0;
    const activeEmbeddingLabels = getActiveEmbeddingLabels(searchTokens, embeddingLabels);
    const displayName = primaryTrace.name === '__count' ? '' : primaryTrace.name;
    return (
        <Box bgcolor={"inherit"} color="inherit" style={{position: 'relative'}}>
            <Box data-testid="chart-extra" color="text.primary" sx={{
                position: 'absolute',
                textAlign: 'right',
                overflow: 'hidden',
                // whiteSpace: 'nowrap',
                textOverflow: 'ellipsis',
                maxWidth: 300,
                right: 8,
                top: 0,
                zIndex: 1000
            }}>
                {displayName !== '' &&
                    <Tooltip title={"Embedding: " + primaryTrace.embedding.name}>
                        <div onClick={handleToggleLegend}
                             style={{cursor: 'pointer', marginRight: 14}}>
                            <Typography color="textPrimary" variant={"subtitle1"}
                                        component={"div"} style={{display: 'inline-block'}}>{displayName}</Typography>
                            <div style={{display: 'inline-block', verticalAlign: 'bottom'}}><ExpandMoreIcon
                                fontSize={"medium"}
                                style={{transform: showLegend ? 'rotate(180deg)' : ''}}/>
                            </div>
                        </div>
                    </Tooltip>
                }

                {primaryTrace.continuous ?
                    <ColorSchemeLegendWrapper
                        handleDomain={onDomain}
                        style={{
                            display: showLegend ? 'block' : 'none'
                        }}
                        handleUpdate={onMeasureFilterUpdated}
                        datasetFilter={datasetFilter}
                        featureSummary={featureSummary}
                        globalFeatureSummary={globalFeatureSummary}
                        nObs={shape[0]}
                        nObsSelected={nObsSelected}
                        name={primaryTrace.name}
                        type={activeFeature.type}
                    /> :
                    <CategoricalLegend
                        legendScrollPosition={legendScrollPosition}
                        handleScrollPosition={handleScrollPosition}
                        visible={showLegend}
                        height={primaryChartSize.height - 40}
                        features={dataset.features}
                        datasetFilter={datasetFilter}
                        handleClick={onDimensionFilterUpdated}
                        handleColorChange={onColorChange}
                        handleNameChange={onCategoricalNameChange}
                        onAddFeatures={onAddFeatures}
                        categoricalNames={categoricalNames}
                        name={primaryTrace.name}
                        scale={primaryTrace.colorScale}
                        nObs={shape[0]}
                        nObsSelected={nObsSelected}
                        globalFeatureSummary={globalFeatureSummary}
                        featureSummary={featureSummary}
                        serverInfo={serverInfo}/>
                }
            </Box>

            {primaryTrace.type === 'scatter' && primaryTrace.embedding.mode == null &&
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
                                   onCamera={onCamera}
                                   handleClick={onDimensionFilterUpdated}

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
                markerOpacity={markerOpacity}
                unselectedMarkerOpacity={unselectedMarkerOpacity}
                onSelected={onSelect}
                onGallery={onGallery}
                handleClick={onDimensionFilterUpdated}
            />}
        </Box>);


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
        legendScrollPosition: state.legendScrollPosition,
        markerOpacity: state.markerOpacity,
        selection: state.selection,
        pointSize: state.pointSize,
        primaryChartSize: state.panel.primaryChartSize,
        searchTokens: state.searchTokens,
        serverInfo: state.serverInfo,
        shape: state.dataset.shape,
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
        handleScrollPosition: value => {
            dispatch(setLegendScrollPosition(value));
        },
        handleWindowSize: () => {
            dispatch(setWindowSize());
        },
        handleEmbeddingData: (value) => {
            dispatch(setEmbeddingData(value));
        },
        handleSearchTokens: (value) => {
            dispatch(setSearchTokens(value, false));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(EmbeddingChart));


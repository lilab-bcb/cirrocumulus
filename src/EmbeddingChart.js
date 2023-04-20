import {Tooltip, Typography} from '@mui/material';
import Box from '@mui/material/Box';
import {find} from 'lodash';
import React, {useEffect, useRef, useState} from 'react';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import PercentIcon from '@mui/icons-material/Percent';
import {connect} from 'react-redux';
import NumbersIcon from '@mui/icons-material/Numbers';
import CheckBoxIcon from '@mui/icons-material/CheckBox';
import {
  getTraceKey,
  handleBrushFilterUpdated,
  handleCategoricalNameChange,
  handleColorChange,
  handleDimensionFilterUpdated,
  handleDomainChange,
  handleMeasureFilterUpdated,
  setCategoricalSortOrder,
  setChartOptions,
  setEmbeddingData,
  setLegendScrollPosition,
  setSearchTokens,
  setWindowSize,
} from './actions';
import CategoricalLegend from './CategoricalLegend';
import ColorSchemeLegendWrapper from './ColorSchemeLegendWrapper';
import ImageChart from './ImageChart';
import MetaEmbedding from './MetaEmbedding';
import ScatterChartThree from './ScatterChartThree';
import {
  FEATURE_TYPE,
  getCategoryValue,
  sortCategories,
  TRACE_TYPE_META_IMAGE,
} from './util';
import memoize from 'memoize-one';
import SortByAlphaIcon from '@mui/icons-material/SortByAlpha';
import IconButton from '@mui/material/IconButton';
import SaveIcon from '@mui/icons-material/Save';
import MenuItem from '@mui/material/MenuItem';
import Menu from '@mui/material/Menu';

const getActiveEmbeddingLabels = memoize((searchTokens, embeddingLabels) => {
  return searchTokens
    .filter(
      (item) =>
        item.type === FEATURE_TYPE.OBS_CAT &&
        embeddingLabels.indexOf(item.id) !== -1,
    )
    .map((item) => item.id);
});

function EmbeddingChart(props) {
  const {
    activeFeature,
    cachedData,
    categoricalNames,
    categoricalSortOrder,
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
    jobResults,
    legendScrollPosition,
    markerOpacity,
    onCategoricalNameChange,
    onChartOptions,
    onColorChange,
    onDimensionFilterUpdated,
    onCategoricalSortOrder,
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
    unselectedPointSize,
  } = props;

  const [showLegend, setShowLegend] = useState(true);
  const [sortOrder, setSortOrder] = useState('alpha'); // alpha, percent, size
  const [selectionMenuOpen, setSelectionMenuOpen] = useState(false);
  const active = 'cirro-active';
  const selectionMenuAnchorEl = useRef();
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

  function handleSortOrder(value) {
    const categories =
      value !== 'alpha'
        ? sortCategories(
            globalFeatureSummary,
            featureSummary,
            activeFeature.name,
            value,
          )
        : null;
    onCategoricalSortOrder({name: activeFeature.name, value: categories});
    setSortOrder(value);
  }

  function handleInvertSelection() {
    setSelectionMenuOpen(false);
    onDimensionFilterUpdated({name: primaryTrace.name, invert: true});
  }

  function handleClearSelection() {
    setSelectionMenuOpen(false);
    onDimensionFilterUpdated({name: primaryTrace.name, clear: true});
  }

  function handleSaveColors() {
    const name = primaryTrace.name;
    const scale = primaryTrace.colorScale;
    const globalDimensionSummary = globalFeatureSummary[name];
    const categories = globalDimensionSummary.categories;
    const renamedCategories = categoricalNames[name] || {};
    const text = [];
    categories.forEach((category) => {
      const renamedCategory = getCategoryValue(renamedCategories, category);
      text.push(renamedCategory);
      text.push('\t');
      text.push(scale(category));
      text.push('\n');
    });
    const blob = new Blob([text.join('')], {
      type: 'text/plain;charset=utf-8',
    });
    window.saveAs(blob, name + '-colors.tsv');
  }

  function onAddFeatures(features) {
    const values = searchTokens.filter(
      (token) => token.type !== FEATURE_TYPE.X,
    );
    const xTokens = searchTokens.filter(
      (token) => token.type === FEATURE_TYPE.X,
    );
    const existingXFeatures = new Set();
    xTokens.forEach((item) => existingXFeatures.add(item.id));
    features.forEach((feature) => {
      if (!existingXFeatures.has(feature)) {
        xTokens.push({id: feature, type: FEATURE_TYPE.X});
        existingXFeatures.add(feature);
      }
    });
    handleSearchTokens(values.concat(xTokens));
  }

  function onCamera(eventName, cameraDef) {
    const primaryTrace = find(
      embeddingData,
      (item) => getTraceKey(item) === activeFeature.embeddingKey,
    );
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
  const primaryTrace = find(
    embeddingData,
    (item) => getTraceKey(item) === activeFeature.embeddingKey,
  );
  if (primaryTrace == null) {
    console.log(activeFeature.embeddingKey + ' not found');
    return null;
  }
  let invertDisabled = true;
  if (!primaryTrace.continuous) {
    let categoricalFilter = datasetFilter[primaryTrace.name];
    if (categoricalFilter != null && categoricalFilter.value.length > 0) {
      let categories;
      for (let i = 0; i < embeddingData.length; i++) {
        if (embeddingData[i].name === primaryTrace.name) {
          categories = embeddingData[i].colorScale.domain();
          break;
        }
      }
      invertDisabled = categories.length === categoricalFilter.value.length;
    }
  }
  const nObsSelected = selection != null ? selection.size : 0;
  const activeEmbeddingLabels = getActiveEmbeddingLabels(
    searchTokens,
    embeddingLabels,
  );
  let displayName;
  if (primaryTrace.featureType === FEATURE_TYPE.JOB_RESULT) {
    let val = find(
      jobResults,
      (jobResult) => jobResult.id === primaryTrace.name,
    );
    displayName = val ? val.name : '';
  } else {
    displayName = primaryTrace.name === '__count' ? '' : primaryTrace.name;
  }
  return (
    <Box bgcolor={'inherit'} color="inherit" style={{position: 'relative'}}>
      <Box
        data-testid="chart-extra"
        color="text.primary"
        sx={{
          position: 'absolute',
          textAlign: 'right',
          right: 8,
          top: 0,
          zIndex: 1000,
        }}
      >
        {displayName !== '' && (
          <>
            {!primaryTrace.continuous && (
              <div style={{float: 'left'}}>
                <Tooltip title={'Sort Legend Alphabetically'}>
                  <span>
                    <IconButton
                      edge={false}
                      disabled={sortOrder !== 'percent'}
                      style={{padding: 0}}
                      className={sortOrder === 'alpha' ? active : ''}
                      size={'small'}
                      aria-label="Sort Alphabetically"
                      onClick={() => handleSortOrder('alpha')}
                    >
                      <SortByAlphaIcon fontSize={'small'} />
                    </IconButton>
                  </span>
                </Tooltip>
                <Tooltip title={'Sort Legend By Percent Selected'}>
                  <IconButton
                    edge={false}
                    size={'small'}
                    style={{padding: 0}}
                    className={sortOrder === 'percent' ? active : ''}
                    aria-label="Sort Legend By Percent Selected"
                    onClick={() => handleSortOrder('percent')}
                  >
                    <PercentIcon />
                  </IconButton>
                </Tooltip>
                <Tooltip title={'Sort Legend By Size'}>
                  <IconButton
                    edge={false}
                    size={'small'}
                    style={{padding: 0}}
                    className={sortOrder === 'size' ? active : ''}
                    aria-label="Sort Legend By Size"
                    onClick={() => handleSortOrder('size')}
                  >
                    <NumbersIcon />
                  </IconButton>
                </Tooltip>

                <Tooltip title={'Selection'}>
                  <span>
                    <IconButton
                      edge={false}
                      size={'small'}
                      ref={selectionMenuAnchorEl}
                      disabled={invertDisabled}
                      style={{padding: 0}}
                      aria-label="Invert Selection"
                      onClick={() => setSelectionMenuOpen(true)}
                    >
                      <CheckBoxIcon />
                    </IconButton>
                    <Menu
                      id="selection-menu"
                      anchorEl={selectionMenuAnchorEl.current}
                      anchorOrigin={{
                        vertical: 'top',
                        horizontal: 'right',
                      }}
                      transformOrigin={{
                        vertical: 'top',
                        horizontal: 'right',
                      }}
                      open={selectionMenuOpen}
                      onClose={() => setSelectionMenuOpen(false)}
                    >
                      <MenuItem onClick={handleClearSelection}>
                        Clear Selection
                      </MenuItem>
                      <MenuItem onClick={handleInvertSelection}>
                        Invert Selection
                      </MenuItem>
                    </Menu>
                  </span>
                </Tooltip>

                <Tooltip title={'Save Colors'}>
                  <IconButton
                    edge={false}
                    size={'small'}
                    style={{
                      paddingLeft: 6,
                      paddingRight: 0,
                      paddingTop: 0,
                      paddingBottom: 0,
                    }}
                    aria-label="Save Colors"
                    onClick={() => handleSaveColors()}
                  >
                    <SaveIcon />
                  </IconButton>
                </Tooltip>
              </div>
            )}
            <Tooltip
              title={
                <>
                  {displayName}
                  <br />
                  {primaryTrace.embedding.name}
                </>
              }
            >
              <div
                onClick={handleToggleLegend}
                style={{
                  cursor: 'pointer',
                  marginRight: 14,
                  display: 'inline-block',
                }}
              >
                <Typography
                  color="textPrimary"
                  variant={'subtitle1'}
                  component={'div'}
                  style={{
                    display: 'inline-block',
                    maxWidth: 200,
                    overflow: 'hidden',
                    whiteSpace: 'nowrap',
                    textOverflow: 'ellipsis',
                  }}
                >
                  {displayName}
                </Typography>
                <div style={{display: 'inline-block', verticalAlign: 'bottom'}}>
                  <ExpandMoreIcon
                    fontSize={'medium'}
                    style={{transform: showLegend ? 'rotate(180deg)' : ''}}
                  />
                </div>
              </div>
            </Tooltip>
          </>
        )}

        {primaryTrace.continuous ? (
          <ColorSchemeLegendWrapper
            handleDomain={onDomain}
            style={{
              display: showLegend ? 'block' : 'none',
            }}
            handleUpdate={onMeasureFilterUpdated}
            datasetFilter={datasetFilter}
            featureSummary={featureSummary}
            globalFeatureSummary={globalFeatureSummary}
            nObs={shape[0]}
            nObsSelected={nObsSelected}
            name={primaryTrace.name}
            type={activeFeature.type}
          />
        ) : (
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
            serverInfo={serverInfo}
            categoricalSortOrder={categoricalSortOrder}
          />
        )}
      </Box>

      {primaryTrace.type === 'scatter' &&
        primaryTrace.embedding.mode == null && (
          <ScatterChartThree
            trace={primaryTrace}
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
          />
        )}

      {primaryTrace.type === TRACE_TYPE_META_IMAGE && (
        <MetaEmbedding
          trace={primaryTrace}
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
        />
      )}
      {primaryTrace.type === 'image' && (
        <ImageChart
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
        />
      )}
    </Box>
  );
}

const mapStateToProps = (state) => {
  return {
    activeFeature: state.activeFeature,
    embeddingData: state.embeddingData,
    cachedData: state.cachedData,
    categoricalNames: state.categoricalNames,
    categoricalSortOrder: state.categoricalSortOrder,
    chartOptions: state.chartOptions,
    dataset: state.dataset,
    datasetFilter: state.datasetFilter,
    embeddingLabels: state.embeddingLabels,
    featureSummary: state.featureSummary,
    globalFeatureSummary: state.globalFeatureSummary,
    legendScrollPosition: state.legendScrollPosition,
    jobResults: state.jobResults,
    markerOpacity: state.markerOpacity,
    selection: state.selection,
    pointSize: state.pointSize,
    primaryChartSize: state.panel.primaryChartSize,
    searchTokens: state.searchTokens,
    serverInfo: state.serverInfo,
    shape: state.dataset.shape,
    unselectedMarkerOpacity: state.unselectedMarkerOpacity,
    unselectedPointSize: state.unselectedPointSize,
  };
};
const mapDispatchToProps = (dispatch) => {
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
    handleScrollPosition: (value) => {
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
    },
    onCategoricalSortOrder: (value) => {
      dispatch(setCategoricalSortOrder(value));
    },
  };
};

export default connect(mapStateToProps, mapDispatchToProps)(EmbeddingChart);

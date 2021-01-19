import {scaleSequential} from 'd3-scale';
import {isPlainObject} from 'lodash';
import {combineReducers} from 'redux';
import {
    ADD_DATASET,
    DELETE_DATASET,
    getTraceKey,
    RESTORE_VIEW,
    SET_BIN_SUMMARY,
    SET_BIN_VALUES,
    SET_CATEGORICAL_COLOR,
    SET_CATEGORICAL_NAME,
    SET_CHART_OPTIONS,
    SET_CHART_SIZE,
    SET_COMBINE_DATASET_FILTERS,
    SET_DATASET,
    SET_DATASET_CHOICES,
    SET_DATASET_FILTER,
    SET_DATASET_FILTERS,
    SET_DIALOG,
    SET_DISTRIBUTION_DATA,
    SET_DISTRIBUTION_PLOT_INTERPOLATOR,
    SET_DISTRIBUTION_PLOT_OPTIONS,
    SET_DOMAIN,
    SET_EMAIL,
    SET_EMBEDDING_DATA,
    SET_EMBEDDING_LABELS,
    SET_FEATURE_SUMMARY,
    SET_GLOBAL_FEATURE_SUMMARY,
    SET_INTERPOLATOR,
    SET_LOADING,
    SET_LOADING_APP,
    SET_MARKER_OPACITY,
    SET_MARKERS,
    SET_MESSAGE,
    SET_NUMBER_OF_BINS,
    SET_POINT_SIZE,
    SET_PRIMARY_CHART_SIZE,
    SET_PRIMARY_TRACE_KEY,
    SET_SAVED_DATASET_STATE,
    SET_SEARCH_TOKENS,
    SET_SELECTED_DISTRIBUTION_DATA,
    SET_SELECTED_EMBEDDING,
    SET_SELECTION,
    SET_SERVER_INFO,
    SET_TAB,
    SET_UNSELECTED_MARKER_OPACITY,
    SET_USER,
    UPDATE_DATASET,
} from '../actions';
import {getInterpolator, updateTraceColors} from '../util';

export const DEFAULT_BIN_SUMMARY = 'max';
export const DEFAULT_NUMBER_BINS = 500;
export const DEFAULT_POINT_SIZE = 1;
export const DEFAULT_MARKER_OPACITY = 1;
export const DEFAULT_UNSELECTED_MARKER_OPACITY = 0.1;
export const DEFAULT_INTERPOLATOR = 'Viridis';
export const DEFAULT_DOT_PLOT_INTERPOLATOR = 'Reds';
export const DEFAULT_DRAG_MODE = 'pan';
export const DEFAULT_SHOW_LABELS = false;
export const DEFAULT_SHOW_AXIS = true;
export const DEFAULT_SHOW_FOG = false;
export const DEFAULT_DARK_MODE = false;
export const DEFAULT_LABEL_FONT_SIZE = 14;
export const DEFAULT_LABEL_STROKE_WIDTH = 4;

const DEFAULT_DIST_PLOT_OPTIONS = {
    chartType: 'dotplot',
    violinScale: 'width',
    violinHeight: 100,
    violinWidth: 80,
    violinShowBoxplot:true
};

const DEFAULT_INTERPOLATOR_OBJ = {
    name: DEFAULT_INTERPOLATOR,
    value: getInterpolator(DEFAULT_INTERPOLATOR)
};

const DEFAULT_DOT_PLOT_INTERPOLATOR_OBJ = {
    name: DEFAULT_DOT_PLOT_INTERPOLATOR,
    value: getInterpolator(DEFAULT_DOT_PLOT_INTERPOLATOR)
};

const DEFAULT_PRIMARY_CHART_SIZE = {
    width: window.innerWidth - 280,
    height: Math.max(300, window.innerHeight - 410)
};

const DEFAULT_CHART_OPTIONS = {
    animating: false,
    dragmode: DEFAULT_DRAG_MODE,
    editSelection: false,
    showAxis: DEFAULT_SHOW_AXIS,
    showFog: DEFAULT_SHOW_FOG,
    darkMode: DEFAULT_DARK_MODE,
    labelFontSize: DEFAULT_LABEL_FONT_SIZE,
    labelStrokeWidth: DEFAULT_LABEL_STROKE_WIDTH
};

function chartSize(state = 300, action) {
    switch (action.type) {
        case SET_CHART_SIZE:
            return action.payload;
        default:
            return state;
    }
}


function primaryChartSize(state = DEFAULT_PRIMARY_CHART_SIZE, action) {
    switch (action.type) {
        case SET_PRIMARY_CHART_SIZE:
            return action.payload;
        default:
            return state;
    }
}

/**
 *
 * @param state Array of value, type where type can be X, obs, or obsCat
 * @param action
 * @returns {*|*[]}
 */
function searchTokens(state = [], action) {
    switch (action.type) {
        case SET_SEARCH_TOKENS:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.q || [];
        default:
            return state;
    }
}

function embeddingLabels(state = [], action) {
    switch (action.type) {
        case SET_EMBEDDING_LABELS:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.embeddingLabels || [];
        default:
            return state;
    }
}

function chartOptions(state = DEFAULT_CHART_OPTIONS, action) {
    switch (action.type) {
        case SET_CHART_OPTIONS:
            return Object.assign({}, action.payload);
        case RESTORE_VIEW:
            return action.payload.chartOptions ? Object.assign(DEFAULT_CHART_OPTIONS, action.payload.chartOptions) : state;
        default:
            return state;
    }
}

function dataset(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
            document.title = action.payload == null ? 'Cirro' : action.payload.name + ' - Cirro';
            return action.payload;
        case UPDATE_DATASET:
            if (action.payload.id === state.id) { // update name, description, url
                document.title = action.payload.name + ' - Cirro';
                return Object.assign(state, action.payload);
            }
        default:
            return state;
    }
}

function datasetChoices(state = [], action) {
    switch (action.type) {
        case SET_DATASET_CHOICES:
            action.payload.sort((a, b) => {
                a = a.name.toLowerCase();
                b = b.name.toLowerCase();
                return a === b ? 0 : (a < b ? -1 : 1);
            });

            return action.payload;
        case SET_EMAIL:
            if (action.payload == null) {
                return [];
            }
            return state;
        case ADD_DATASET:
            state.push(action.payload);
            state.sort((a, b) => {
                a = a.name.toLowerCase();
                b = b.name.toLowerCase();
                return a === b ? 0 : (a < b ? -1 : 1);
            });
            return state.slice();
        case UPDATE_DATASET:
        case DELETE_DATASET:
            let index = -1;
            for (let i = 0; i < state.length; i++) {
                if (state[i].id === action.payload.id) {
                    index = i;
                    break;
                }
            }
            if (index !== -1) {
                if (action.type === UPDATE_DATASET) {
                    state[index] = action.payload;
                    state.sort((a, b) => {
                        a = a.name.toLowerCase();
                        b = b.name.toLowerCase();
                        return a === b ? 0 : (a < b ? -1 : 1);
                    });
                } else {
                    state.splice(index, 1);
                }

                return state.slice();
            }
            return state;
        default:
            return state;
    }
}

// set the selected embeddings, each embedding has name (str) e.g X_umap, nbins (int), _nbins (str), agg (str), bin (boolean), dimensions (int), precomputed (bool)
function embeddings(state = [], action) {
    switch (action.type) {
        case SET_SELECTED_EMBEDDING:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.embeddings != null ? action.payload.embeddings : state;
        default:
            return state;
    }
}

function markers(state = [], action) {
    switch (action.type) {
        case SET_MARKERS:
            return action.payload;
        case SET_DATASET:
            let result = action.payload == null ? [] : action.payload.markers || [];
            if (isPlainObject(result)) { // old style, category=>name=>features
                let newResults = [];
                for (let categoryName in result) {
                    let category = result[categoryName];
                    for (let name in category) {
                        newResults.push({
                            category: categoryName,
                            name: name,
                            id: categoryName + '-' + name,
                            readonly: true,
                            features: category[name]
                        });
                    }
                }
                result = newResults;
            }

            if (action.payload != null && action.payload.markers_read_only != null) {
                if (isPlainObject(action.payload.markers_read_only)) { // old style, name => features
                    let markers = action.payload.markers_read_only;
                    for (let categoryName in markers) {
                        let category = markers[categoryName];
                        for (let name in category) {
                            result.push({
                                category: categoryName,
                                name: name,
                                id: categoryName + '-' + name,
                                readonly: true,
                                features: category[name]
                            });
                        }
                    }
                } else {
                    action.payload.markers_read_only.forEach(item => {
                        item.readonly = true;
                        result.push(item);
                    });
                }
            }
            result.forEach(item => {
                if (item.id == null) {
                    item.id = item.category + '-' + item.name;
                }
            });

            return result;
        default:
            return state;
    }
}


function binSummary(state = DEFAULT_BIN_SUMMARY, action) {
    switch (action.type) {
        case SET_BIN_SUMMARY:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_BIN_SUMMARY;
        case RESTORE_VIEW:
            return action.payload.binSummary != null ? action.payload.binSummary : state;
        default:
            return state;
    }
}


function binValues(state = false, action) {
    switch (action.type) {
        case SET_BIN_VALUES:
            return action.payload;
        case SET_DATASET:
            return false;
        case RESTORE_VIEW:
            return action.payload.binValues != null ? action.payload.binValues : state;
        default:
            return state;
    }
}


function numberOfBins(state = DEFAULT_NUMBER_BINS, action) {
    switch (action.type) {
        case SET_NUMBER_OF_BINS:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_NUMBER_BINS;
        case RESTORE_VIEW:
            return action.payload.numberOfBins != null ? action.payload.numberOfBins : state;
        default:
            return state;
    }
}


function markerOpacity(state = DEFAULT_MARKER_OPACITY, action) {
    switch (action.type) {
        case SET_MARKER_OPACITY:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_MARKER_OPACITY;
        case RESTORE_VIEW:
            return action.payload.markerOpacity != null ? action.payload.markerOpacity : state;
        default:
            return state;
    }
}

function unselectedMarkerOpacity(state = DEFAULT_UNSELECTED_MARKER_OPACITY, action) {
    switch (action.type) {
        case SET_UNSELECTED_MARKER_OPACITY:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_UNSELECTED_MARKER_OPACITY;
        case RESTORE_VIEW:
            return action.payload.unselectedMarkerOpacity != null ? action.payload.unselectedMarkerOpacity : state;
        default:
            return state;
    }
}


function message(state = null, action) {
    switch (action.type) {
        case SET_MESSAGE:
            return action.payload;
        default:
            return state;
    }
}

function email(state = null, action) {
    switch (action.type) {
        case SET_EMAIL:
            return action.payload;
        default:
            return state;
    }
}

function user(state = {}, action) {
    switch (action.type) {
        case SET_USER:
            return action.payload;
        default:
            return state;
    }
}

/**
 * Object that contains count (number), chart (object). Each key in chart is the full layout name. Each value contains
 * userPoints (selected points in chart space) and points (the selected points in bin space if binning)
 */
function selection(state = {chart: {}}, action) {
    switch (action.type) {
        case SET_SELECTION:
            return action.payload;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

/**
 * Feature summary maps measure and dimension names to another object containing
 * categories and counts for dimensions, and statistics such as min, max for measures.
 * Features summaries are in the space of selected cells.
 */
function featureSummary(state = {}, action) {
    switch (action.type) {
        case SET_FEATURE_SUMMARY:
            return action.payload;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

/**
 * Feature summary maps measure and dimension names to another object containing
 * categories and counts for dimensions, and statistics such as min, max for measures.
 * Features summaries are in the space of all cells.
 */
function globalFeatureSummary(state = {}, action) {
    switch (action.type) {
        case SET_GLOBAL_FEATURE_SUMMARY:
            return action.payload;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

function serverInfo(state = {}, action) {
    switch (action.type) {
        case SET_SERVER_INFO:
            return action.payload;
        default:
            return state;
    }
}


function loadingApp(state = {loading: false, progress: 0}, action) {
    switch (action.type) {
        case SET_LOADING_APP:
            return action.payload;
        default:
            return state;
    }
}


function tab(state = 'embedding', action) {
    switch (action.type) {
        case SET_TAB:
            return action.payload;
        default:
            return state;
    }
}

function dialog(state = null, action) {
    switch (action.type) {
        case SET_DIALOG:
            return action.payload;
        default:
            return state;
    }
}


function dotPlotInterpolator(state = DEFAULT_DOT_PLOT_INTERPOLATOR_OBJ, action) {
    switch (action.type) {
        case SET_DISTRIBUTION_PLOT_INTERPOLATOR:
            return action.payload;
        case RESTORE_VIEW:
            return action.payload.dotPlotInterpolator != null ? action.payload.dotPlotInterpolator : state;
        default:
            return state;
    }
}


/*
sortBy, minSize, maxSize, min, max, chartType
 */
function distributionPlotOptions(state = DEFAULT_DIST_PLOT_OPTIONS, action) {
    switch (action.type) {
        case SET_DATASET:
            return DEFAULT_DIST_PLOT_OPTIONS;
        case SET_DISTRIBUTION_PLOT_OPTIONS:
            return Object.assign({}, state, action.payload);
        case RESTORE_VIEW:
            return action.payload.distributionPlotOptions != null ? action.payload.distributionPlotOptions : state;
        default:
            return state;
    }
}

function distributionData(state = [], action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            return state.slice();
        case SET_DISTRIBUTION_DATA:
            return action.payload;
        case SET_DATASET:
            return [];
        default:
            return state;
    }
}

function selectedDistributionData(state = [], action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            return state.slice();
        case SET_SELECTED_DISTRIBUTION_DATA:
            return action.payload;
        case SET_DATASET:
            return [];
        default:
            return state;
    }
}


function datasetFilters(state = [], action) {
    switch (action.type) {
        case SET_DATASET:
            return [];
        case SET_DATASET_FILTERS:
            //action.payload.sort()
            return action.payload;
        default:
            return state;
    }
}


function datasetFilter(state = {}, action) {
    switch (action.type) {
        case SET_DATASET:
            return {};
        case SET_DATASET_FILTER:
            return action.payload;
        case RESTORE_VIEW:
            return action.payload.datasetFilter != null ? action.payload.datasetFilter : state;
        default:
            return state;
    }
}

function combineDatasetFilters(state = 'and', action) {
    switch (action.type) {
        case SET_COMBINE_DATASET_FILTERS:
            return action.payload;
        case SET_DATASET:
            return "and";
        case RESTORE_VIEW:
            return action.payload.combineDatasetFilters != null ? action.payload.combineDatasetFilters : state;
        default:
            return state;
    }
}

// category -> value -> newValue
function categoricalNames(state = {}, action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            let category = state[action.payload.name];
            if (category === undefined) {
                category = {};
                state[action.payload.name] = category;
            }
            if (action.payload.value == null || action.payload.value === '') {
                delete category[action.payload.oldValue];
            } else {
                category[action.payload.oldValue] = action.payload.value;
            }
            return Object.assign({}, state);
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

// maps keys to values, features=>array, embedding key=>x,y,z, feature_embedding_key=>array for binned embeddings
function cachedData(state = {}, action) {
    switch (action.type) {
        case SET_DATASET :
            return {};
    }
    return state;

}

// each item has  data (list of traces, each trace has x, y, etc.), layout
function embeddingData(state = [], action) {
    switch (action.type) {
        case SET_EMBEDDING_DATA :
            return action.payload;
        case SET_SELECTION:
            return state.slice();
        case SET_DOMAIN:
            state.forEach((traceInfo, stateIndex) => {
                if (traceInfo.continuous && traceInfo.name === action.payload.name) {
                    const summary = action.payload.summary;
                    let domain = [summary.min, summary.max];
                    if (summary.customMin != null && !isNaN(summary.customMin)) {
                        domain[0] = summary.customMin;
                    }
                    if (summary.customMax != null && !isNaN(summary.customMax)) {
                        domain[1] = summary.customMax;
                    }
                    traceInfo.colorScale.domain(domain);
                    updateTraceColors(traceInfo);
                }
            });
            return state.slice();
        case SET_CATEGORICAL_COLOR:
            state.forEach((traceInfo, stateIndex) => {
                if (!traceInfo.continuous && traceInfo.name === action.payload.name) {
                    let index = traceInfo.colorScale.domain().indexOf(action.payload.value);
                    let range = traceInfo.colorScale.range();
                    range[index] = action.payload.color;
                    traceInfo.colorScale.range(range);
                    updateTraceColors(traceInfo);
                }
            });
            return state.slice();

        case SET_INTERPOLATOR:
            // update colors for existing continuous traces

            state.forEach((traceInfo, stateIndex) => {
                if (traceInfo.continuous) {
                    let domain = traceInfo.colorScale.domain();
                    traceInfo.colorScale = scaleSequential(action.payload.value).domain(domain);
                    updateTraceColors(traceInfo);
                }
            });
            return state.slice();
        case SET_DATASET:
            return [];
        default:
            return state;
    }
}

export function pointSize(state = DEFAULT_POINT_SIZE, action) {
    switch (action.type) {
        case SET_POINT_SIZE:
            return action.payload;
        case RESTORE_VIEW:
            return action.payload.pointSize != null ? action.payload.pointSize : state;
        default:
            return state;
    }
}


// used to restore state when toggling datasets
export function savedDatasetState(state = {}, action) {
    switch (action.type) {
        case SET_SAVED_DATASET_STATE:
            return action.payload;
        default:
            return state;
    }
}

export function primaryTraceKey(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
            return null;
        case SET_PRIMARY_TRACE_KEY:
            return action.payload;
        case SET_EMBEDDING_DATA:
            let traces = action.payload.filter(traceInfo => traceInfo.active);
            if (traces.length === 0) {
                return null;
            }
            return getTraceKey(traces[traces.length - 1]); // last feature becomes primary
        default:
            return state;
    }
}

function loading(state = false, action) {
    switch (action.type) {
        case SET_LOADING:
            return action.payload;
        default:
            return state;
    }
}


function interpolator(state = DEFAULT_INTERPOLATOR_OBJ, action) {
    switch (action.type) {
        case SET_INTERPOLATOR:
            return action.payload;
        case RESTORE_VIEW:
            return action.payload.colorScheme != null ? action.payload.colorScheme : state;
        default:
            return state;
    }
}

export default combineReducers({
    binSummary,
    binValues,
    cachedData,
    categoricalNames,
    chartOptions,
    chartSize,
    combineDatasetFilters,
    dataset,
    datasetChoices,
    datasetFilter,
    datasetFilters,
    dialog,
    distributionData,
    distributionPlotOptions,
    dotPlotInterpolator,
    email,
    embeddingData,
    embeddingLabels,
    embeddings,
    featureSummary,
    globalFeatureSummary,
    interpolator,
    loading,
    loadingApp,
    markerOpacity,
    markers,
    message,
    numberOfBins,
    pointSize,
    primaryChartSize,
    primaryTraceKey,
    savedDatasetState,
    searchTokens,
    selectedDistributionData,
    selection,
    serverInfo,
    tab,
    unselectedMarkerOpacity,
    user

});

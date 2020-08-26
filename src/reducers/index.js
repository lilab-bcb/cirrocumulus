import {scaleSequential} from 'd3-scale';
import {combineReducers} from 'redux';
import {
    ADD_DATASET,
    DELETE_DATASET,
    DIFF_EXP_RESULTS,
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
    SET_DOT_PLOT_DATA,
    SET_DOT_PLOT_SORT_ORDER,
    SET_EMAIL,
    SET_EMBEDDING_DATA,
    SET_FEATURE_SUMMARY,
    SET_GLOBAL_FEATURE_SUMMARY,
    SET_INTERPOLATOR,
    SET_LOADING,
    SET_LOADING_APP,
    SET_MARKER_OPACITY,
    SET_MARKER_OPACITY_UI,
    SET_MESSAGE,
    SET_NUMBER_OF_BINS,
    SET_NUMBER_OF_BINS_UI,
    SET_POINT_SIZE,
    SET_PRIMARY_TRACE_KEY,
    SET_SEARCH_TOKENS,
    SET_SELECTED_DOT_PLOT_DATA,
    SET_SELECTED_EMBEDDING,
    SET_SELECTION,
    SET_SERVER_INFO,
    SET_TAB,
    SET_UNSELECTED_MARKER_OPACITY,
    SET_UNSELECTED_MARKER_OPACITY_UI,
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
export const DEFAULT_DRAG_MODE = 'pan';
export const DEFAULT_SHOW_LABELS = false;
export const DEFAULT_SHOW_AXIS = true;
export const DEFAULT_SHOW_FOG = true;
export const DEFAULT_DARK_MODE = false;

const DEFAULT_INTERPOLATOR_OBJ = {
    name: DEFAULT_INTERPOLATOR,
    reversed: false,
    value: getInterpolator(DEFAULT_INTERPOLATOR)
};

const DEFAULT_CHART_OPTIONS = {
    animating: false, dragmode: DEFAULT_DRAG_MODE, editSelection: false, showLabels: DEFAULT_SHOW_LABELS,
    showAxis: DEFAULT_SHOW_AXIS, showFog: DEFAULT_SHOW_FOG, darkMode: DEFAULT_DARK_MODE
};

function chartSize(state = 500, action) {
    switch (action.type) {
        case SET_CHART_SIZE:
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

function diffExpResults(state = null, action) {
    switch (action.type) {
        case DIFF_EXP_RESULTS:
            return action.payload;
        case SET_DATASET:
            return null;
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

function numberOfBinsUI(state = DEFAULT_NUMBER_BINS, action) {
    switch (action.type) {
        case SET_NUMBER_OF_BINS:
        case SET_NUMBER_OF_BINS_UI:
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

function markerOpacityUI(state = 1, action) {
    switch (action.type) {
        case SET_MARKER_OPACITY:
        case SET_MARKER_OPACITY_UI:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_MARKER_OPACITY;
        case RESTORE_VIEW:
            return action.payload.markerOpacity != null ? action.payload.markerOpacity : state;
        default:
            return state;
    }
}

function unselectedMarkerOpacityUI(state = DEFAULT_UNSELECTED_MARKER_OPACITY, action) {
    switch (action.type) {
        case SET_UNSELECTED_MARKER_OPACITY:
        case SET_UNSELECTED_MARKER_OPACITY_UI:
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

function datasetChoices(state = [], action) {
    switch (action.type) {
        case SET_DATASET_CHOICES:
            return action.payload;
        case SET_EMAIL:
            if (action.payload == null) {
                return [];
            }
            return state;
        case ADD_DATASET:
            state = state.concat([action.payload]);
            state.sort((a, b) => a.name.toLowerCase() < b.name.toLowerCase() ? -1 : 1);
            return state;
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

// [{
//     "categories": array,
//     "name": str,
//     "fractionRange":[number],
//     "meanRange":[number]
//     "values": [{
//         "name": str,
//         "fractionExpressed": array
//         "mean": array,
//         "active": bool
//     }]
// }]

function updateDotPlotDataRange(data) {
    data.forEach(categoryItem => {
        let fractionRange = [Number.MAX_VALUE, -Number.MAX_VALUE];
        let meanRange = [Number.MAX_VALUE, -Number.MAX_VALUE];
        categoryItem.values.forEach(feature => {
            for (let i = 0, n = feature.mean.length; i < n; i++) {
                fractionRange[0] = Math.min(feature.fractionExpressed[i], fractionRange[0]);
                fractionRange[1] = Math.max(feature.fractionExpressed[i], fractionRange[1]);
                meanRange[0] = Math.min(feature.mean[i], meanRange[0]);
                meanRange[1] = Math.max(feature.mean[i], meanRange[1]);
            }
        });
        categoryItem.meanRange = meanRange;
        categoryItem.fractionRange = fractionRange;
    });
}

function dotPlotData(state = [], action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            for (let i = 0; i < state.length; i++) {
                let item = state[i];
                state[i] = Object.assign({}, item);
            }
            return state.slice();
        case SET_DOT_PLOT_SORT_ORDER:
            const name = action.payload.name;
            const sortBy = action.payload.value;
            for (let i = 0; i < state.length; i++) {
                let item = state[i];
                if (item.name === name) {
                    item.sortBy = sortBy;
                    state[i] = Object.assign({}, item);
                    break;
                }
            }
            return state.slice();
        case SET_DOT_PLOT_DATA:
            updateDotPlotDataRange(action.payload);
            return action.payload;
        case SET_DATASET:
            return [];
        default:
            return state;
    }
}

function selectedDotPlotData(state = [], action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            for (let i = 0; i < state.length; i++) {
                let item = state[i];
                state[i] = Object.assign({}, item);
            }
            return state.slice();
        case SET_DOT_PLOT_SORT_ORDER:
            const name = action.payload.name;
            const sortBy = action.payload.value;
            for (let i = 0; i < state.length; i++) {
                let item = state[i];
                if (item.name === name) {
                    item.sortBy = sortBy;
                    state[i] = Object.assign({}, item);
                    break;
                }
            }
            return state.slice();
        case SET_SELECTED_DOT_PLOT_DATA:
            updateDotPlotDataRange(action.payload);
            return action.payload;
        case SET_DATASET:
            return [];
        default:
            return state;
    }
}


function datasetFilters(state = [], action) {
    switch (action.type) {
        case SET_DATASET_FILTERS:
            //action.payload.sort()
            return action.payload;
        default:
            return state;
    }
}


function datasetFilter(state = {}, action) {
    switch (action.type) {
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


// each item has  data (list of traces, each trace has x, y, etc.), layout
function embeddingData(state = [], action) {
    switch (action.type) {
        case SET_EMBEDDING_DATA :
            return action.payload;
        case SET_SELECTION:
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
                    // ensure domain is in order min to max if not reversed, otherwise max to min
                    if (action.payload.reversed) {
                        domain.sort((a, b) => b - a);
                    } else {
                        domain.sort((a, b) => a - b);
                    }

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
        default:
            return state;
    }
}

export function primaryTraceKey(state = null, action) {
    switch (action.type) {
        case SET_PRIMARY_TRACE_KEY:
            return action.payload;
        case SET_EMBEDDING_DATA:
            let traces = action.payload.filter(traceInfo => traceInfo.active);
            if (traces.length === 0) {
                return null;
            }
            const activeTrace = traces[traces.length - 1];
            return getTraceKey(activeTrace); // last feature becomes primary
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
    categoricalNames,
    chartOptions,
    chartSize,
    combineDatasetFilters,
    datasetFilter,
    datasetFilters,
    dataset,
    datasetChoices,
    dialog,
    dotPlotData,
    diffExpResults,
    email,
    embeddingData,
    embeddings,
    featureSummary,
    globalFeatureSummary,
    interpolator,
    loading,
    loadingApp,
    markerOpacity,
    markerOpacityUI,
    numberOfBins,
    numberOfBinsUI,
    message,
    pointSize,
    primaryTraceKey,
    searchTokens,
    selection,
    selectedDotPlotData,
    serverInfo,
    tab,
    unselectedMarkerOpacity,
    unselectedMarkerOpacityUI,
    user
});

import {color} from 'd3-color';
import {scaleSequential} from 'd3-scale';
import {combineReducers} from 'redux';
import {
    ADD_DATASET,
    DELETE_DATASET,
    getTraceKey,
    RESTORE_VIEW,
    SET_BIN_SUMMARY,
    SET_BIN_VALUES,
    SET_CATEGORICAL_COLOR,
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
    SET_FEATURES,
    SET_FEATURES_UI,
    SET_GLOBAL_FEATURE_SUMMARY,
    SET_GROUP_BY,
    SET_INTERPOLATOR,
    SET_LOADING,
    SET_LOADING_APP,
    SET_MARKER_OPACITY,
    SET_MARKER_OPACITY_UI,
    SET_MESSAGE,
    SET_NUMBER_OF_BINS,
    SET_NUMBER_OF_BINS_UI,
    SET_PRIMARY_TRACE_KEY,
    SET_SELECTED_EMBEDDING,
    SET_SELECTION,
    SET_SERVER_INFO,
    SET_UNSELECTED_MARKER_OPACITY,
    SET_UNSELECTED_MARKER_OPACITY_UI,
    SET_USER,
    UPDATE_DATASET,
} from '../actions';
import {getInterpolator, getRgbScale} from '../PlotUtil';

export const DEFAULT_BIN_SUMMARY = 'max';
export const DEFAULT_NUMBER_BINS = 500;

export const DEFAULT_MARKER_OPACITY = 1;
export const DEFAULT_UNSELECTED_MARKER_OPACITY = 0.1;
export const DEFAULT_INTERPOLATOR = 'Viridis';
const DEFAULT_INTERPOLATOR_OBJ = {
    name: DEFAULT_INTERPOLATOR,
    reversed: false,
    value: getInterpolator(DEFAULT_INTERPOLATOR)
};

function features(state = [], action) {
    switch (action.type) {
        case SET_FEATURES:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.features || [];
        default:
            return state;
    }
}

function featuresUI(state = [], action) {
    switch (action.type) {
        case SET_FEATURES:
        case SET_FEATURES_UI:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.features || [];
        default:
            return state;
    }
}


function groupBy(state = [], action) {
    switch (action.type) {
        case SET_GROUP_BY:
            return action.payload;
        case SET_DATASET:
            return [];
        case RESTORE_VIEW:
            return action.payload.groupBy || [];
        default:
            return state;
    }
}


function dataset(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
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

function dialog(state = null, action) {
    switch (action.type) {
        case SET_DIALOG:
            return action.payload;
        default:
            return state;
    }
}

function dotPlotData(state = [], action) {
    switch (action.type) {
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

function updateChartColorScale(traceInfo) {
    const rgbScale = getRgbScale();
    let colorScale = traceInfo.colorScale;
    let colors = [];
    const colorMapper = traceInfo.isImage ? rgb => rgb.formatHex() : rgb => {
        return {
            r: rgbScale(rgb.r),
            g: rgbScale(rgb.g),
            b: rgbScale(rgb.b),
            opacity: 1
        };
    };
    for (let i = 0, n = traceInfo.x.length; i < n; i++) {
        let rgb = color(colorScale(traceInfo.values[i]));
        colors.push(colorMapper(rgb));
    }
    traceInfo.marker.color = colors;
    traceInfo.colorScale = colorScale;
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
                    updateChartColorScale(traceInfo);
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
                    updateChartColorScale(traceInfo);
                }

            });
            return state.slice();
        case SET_DATASET:
            return [];
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
            return getTraceKey(activeTrace);
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
    combineDatasetFilters,
    datasetFilter,
    datasetFilters,
    dataset,
    datasetChoices,
    dialog,
    dotPlotData,
    email,
    embeddingData,
    embeddings,
    features,
    featureSummary,
    featuresUI,
    globalFeatureSummary,
    groupBy,
    interpolator,
    loading,
    loadingApp,
    markerOpacity,
    markerOpacityUI,
    numberOfBins,
    numberOfBinsUI,
    message,
    primaryTraceKey,
    selection,
    serverInfo,
    unselectedMarkerOpacity,
    unselectedMarkerOpacityUI,
    user
});

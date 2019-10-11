import {color} from 'd3-color';
import {scaleLinear, scaleSequential} from 'd3-scale';
import {combineReducers} from 'redux';
import {
    ADD_DATASET,
    DELETE_DATASET,
    RESTORE_VIEW,
    SET_BIN_SUMMARY,
    SET_BIN_VALUES,
    SET_DATASET,
    SET_DATASET_CHOICES,
    SET_DIALOG,
    SET_DOT_PLOT_DATA,
    SET_EMAIL,
    SET_EMBEDDING_DATA,
    SET_FEATURES,
    SET_GROUP_BY,
    SET_INTERPOLATOR,
    SET_LOADING,
    SET_LOADING_APP,
    SET_MARKER_OPACITY,
    SET_MARKER_OPACITY_UI,
    SET_MARKER_SIZE,
    SET_MARKER_SIZE_UI,
    SET_MESSAGE,
    SET_NUMBER_OF_BINS,
    SET_NUMBER_OF_BINS_UI,
    SET_SELECTED_EMBEDDING,
    SET_SELECTED_POINTS,
    SET_SERVER_INFO,
    SET_USER,
    UPDATE_DATASET,
} from '../actions';
import PlotUtil, {getInterpolator} from '../PlotUtil';

export const DEFAULT_MARKER_SIZE = 5;
export const DEFAULT_MARKER_OPACITY = 1;
export const DEFAULT_BIN_SUMMARY = 'mean';
export const DEFAULT_NUMBER_BINS = 500;
export const DEFAULT_INTERPOLATOR = 'Viridis';
const DEFAULT_INTERPOLATOR_OBJ = {name: DEFAULT_INTERPOLATOR, value: getInterpolator(DEFAULT_INTERPOLATOR)};

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

// set the selected embeddings
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

function dataset(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
            return action.payload;
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

function markerSize(state = DEFAULT_MARKER_SIZE, action) {
    switch (action.type) {
        case SET_MARKER_SIZE:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_MARKER_SIZE;
        case RESTORE_VIEW:
            return action.payload.markerSize != null ? action.payload.markerSize : state;
        default:
            return state;
    }
}

function selectedpoints(state = [], action) {
    switch (action.type) {
        case SET_SELECTED_POINTS:
            return action.payload;
        case SET_DATASET:
            return [];
        case SET_SELECTED_EMBEDDING:
            return []; // TODO persist
        case RESTORE_VIEW:
            return action.payload.selectedpoints != null ? action.payload.selectedpoints : state;
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

function markerSizeUI(state = null, action) {
    switch (action.type) {
        case SET_MARKER_SIZE:
        case SET_MARKER_SIZE_UI:
            return action.payload;
        case SET_DATASET:
            return DEFAULT_MARKER_SIZE;
        case RESTORE_VIEW:
            return action.payload.markerSize != null ? action.payload.markerSize : state;
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


function serverInfo(state = {}, action) {
    switch (action.type) {
        case SET_SERVER_INFO:
            return action.payload;
        default:
            return state;
    }
}


function loadingApp(state = true, action) {
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

function dotPlotData(state = null, action) {
    switch (action.type) {
        case SET_DOT_PLOT_DATA:
            return action.payload;
        case SET_DATASET:
            return null;
        default:
            return state;
    }
}

function embeddingData(state = [], action) {
    switch (action.type) {
        case SET_EMBEDDING_DATA :
            return action.payload;
        case SET_MARKER_SIZE:
            state.forEach(item => {
                item.data.forEach(trace => {
                    trace.marker.size = action.payload;
                });
                item.data = item.data.slice();
            });
            return state.slice();
        case SET_SELECTED_POINTS:
            state.forEach(item => {
                item.data.forEach(trace => {
                    trace.selectedpoints = action.payload;
                });

                item.data = item.data.slice();
            });
            return state.slice();
        case SET_INTERPOLATOR:
            // update colors for existing traces
            // TODO custom categorical colors
            let rgbScale = scaleLinear().domain([0, 255]).range([0, 1]);
            state.forEach(item => {
                item.data.forEach(trace => {
                    if (trace.domain) {
                        let colorScale = scaleSequential(action.payload.value).domain(trace.domain);
                        let colors = [];
                        for (let i = 0, n = trace.values.length; i < n; i++) {
                            let rgb = color(colorScale(trace.values[i]));
                            colors.push([rgbScale(rgb.r), rgbScale(rgb.g), rgbScale(rgb.b)]);
                        }
                        trace.marker.color = colors;
                    }
                });
                item.data = item.data.slice();
            });
            return state.slice();
        case SET_MARKER_OPACITY:
            state.forEach(item => {
                item.data.forEach(trace => {
                    trace.marker.opacity = action.payload;
                });
                item.data = item.data.slice();
            });
            return state.slice();
        case SET_DATASET:
            return [];
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

function plotConfig(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
            return PlotUtil.createPlotConfig();
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
    features,
    groupBy,
    embeddings,
    numberOfBins,
    binSummary,
    email,
    user,
    loadingApp,
    dataset,
    binValues,
    datasetChoices,
    dialog,
    dotPlotData,
    embeddingData,
    markerSize,
    loading,
    plotConfig,
    message,
    selectedpoints,
    markerOpacity,
    markerOpacityUI,
    numberOfBinsUI,
    markerSizeUI,
    interpolator,
    serverInfo
});

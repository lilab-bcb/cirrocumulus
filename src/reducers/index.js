import {isPlainObject, isString} from 'lodash';
import {combineReducers} from 'redux';
import {
    ADD_DATASET,
    ADD_TASK,
    DEFAULT_DARK_MODE,
    DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR,
    DEFAULT_DRAG_MODE,
    DEFAULT_INTERPOLATORS,
    DEFAULT_LABEL_FONT_SIZE,
    DEFAULT_LABEL_STROKE_WIDTH,
    DEFAULT_MARKER_OPACITY,
    DEFAULT_POINT_SIZE,
    DEFAULT_SHOW_AXIS,
    DEFAULT_SHOW_FOG,
    DEFAULT_UNSELECTED_MARKER_OPACITY,
    DELETE_DATASET,
    REMOVE_TASK,
    RESTORE_VIEW,
    SET_ACTIVE_FEATURE,
    SET_CATEGORICAL_NAME,
    SET_CHART_OPTIONS,
    SET_CHART_SIZE,
    SET_COMBINE_DATASET_FILTERS,
    SET_DATASET,
    SET_DATASET_CHOICES,
    SET_DATASET_FILTER,
    SET_DATASET_FILTERS,
    SET_DATASET_VIEWS,
    SET_DIALOG,
    SET_DISTRIBUTION_DATA,
    SET_DISTRIBUTION_PLOT_INTERPOLATOR,
    SET_DISTRIBUTION_PLOT_OPTIONS,
    SET_DOMAIN,
    SET_DRAG_DIVIDER,
    SET_DRAWER_OPEN,
    SET_EMAIL,
    SET_EMBEDDING_DATA,
    SET_EMBEDDING_LABELS,
    SET_FEATURE_SUMMARY,
    SET_GLOBAL_FEATURE_SUMMARY,
    SET_INTERPOLATOR,
    SET_JOB_RESULT,
    SET_JOB_RESULTS,
    SET_LEGEND_SCROLL_POSITION,
    SET_LOADING_APP,
    SET_MARKER_OPACITY,
    SET_MARKERS,
    SET_MESSAGE,
    SET_POINT_SIZE,
    SET_SAVED_DATASET_STATE,
    SET_SEARCH_TOKENS,
    SET_SELECTED_DISTRIBUTION_DATA,
    SET_SELECTED_EMBEDDING,
    SET_SELECTION,
    SET_SERVER_INFO,
    SET_TAB,
    SET_UNSELECTED_MARKER_OPACITY,
    SET_UNSELECTED_POINT_SIZE,
    SET_USER,
    SET_WINDOW_SIZE,
    UPDATE_CATEGORICAL_COLOR,
    UPDATE_CATEGORICAL_NAME,
    UPDATE_DATASET
} from '../actions';
import {createCategoryToStats} from '../MetaEmbedding';
import {
    createColorScale,
    FEATURE_TYPE,
    getInterpolator,
    INTERPOLATOR_SCALING_NONE,
    NATSORT,
    TRACE_TYPE_META_IMAGE,
    updateTraceColors
} from '../util';


export const DIST_PLOT_OPTIONS = {
    chartType: 'dotplot',
    violinScale: 'width',
    violinHeight: 100,
    violinWidth: 80,
    violinShowBoxplot: true
};
const DEFAULT_DIST_PLOT_OPTIONS = {
    X: DIST_PLOT_OPTIONS,
    modules: Object.assign({}, DIST_PLOT_OPTIONS, {chartType: 'heatmap'}),
    obs: Object.assign({}, DIST_PLOT_OPTIONS, {chartType: 'violin'})
};


export const DISTRIBUTION_PLOT_INTERPOLATOR_OBJ = {
    name: DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR,
    value: getInterpolator(DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR),
    reversed: false,
    scale: INTERPOLATOR_SCALING_NONE
};

const DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR_OBJ = {
    X: Object.assign({}, DISTRIBUTION_PLOT_INTERPOLATOR_OBJ),
    modules: Object.assign({}, DISTRIBUTION_PLOT_INTERPOLATOR_OBJ),
    obs: Object.assign({}, DISTRIBUTION_PLOT_INTERPOLATOR_OBJ)
};
const DEFAULT_PRIMARY_CHART_SIZE = {
    width: window.innerWidth - 280,
    height: Math.max(300, window.innerHeight - 370)
};

const DEFAULT_CHART_OPTIONS = {
    animating: false,
    dragmode: DEFAULT_DRAG_MODE,
    showGalleryLabels: false,
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


function panel(state = {dividerDelta: 0, drawerOpen: true, primaryChartSize: DEFAULT_PRIMARY_CHART_SIZE}, action) {
    switch (action.type) {
        case SET_DRAG_DIVIDER:
        case SET_DRAWER_OPEN:
        case SET_WINDOW_SIZE:
            let primaryChartSize = state.primaryChartSize;
            const drawerOpen = action.type === SET_DRAWER_OPEN ? action.payload : state.drawerOpen;

            if (action.type === SET_WINDOW_SIZE) {
                const windowWidth = window.innerWidth;
                const windowHeight = window.innerHeight;
                const width = windowWidth - (drawerOpen ? 280 : 40);
                const height = Math.max(300, windowHeight - 370);
                primaryChartSize = {width: width, height: height};
            } else if (action.type === SET_DRAG_DIVIDER) {
                primaryChartSize = {width: primaryChartSize.width, height: action.payload};
            }
            return {drawerOpen: drawerOpen, primaryChartSize: primaryChartSize};
        default:
            return state;
    }
}


/**
 *
 * @param state Array of {id:str, type:str} where type is one of FEATURE_TYPE
 * @param action
 * @returns Array of search token objects
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
            if (action.payload) {
                let features = action.payload.var;
                if (features) {
                    if (features.length > 0 && isString(features[0])) {
                        features = features.map(item => {
                            return {id: item, group: ''};
                        });
                    }
                    features.forEach(item => {
                        if (item.group == null) {
                            item.group = ''; // set default group
                        }
                        item.text = item.id;
                        if (item.text.startsWith(item.group + '-')) { // hide group
                            item.text = item.text.substring(item.group.length + 1);
                        }
                    });
                    features.sort((item1, item2) => {
                        const g = NATSORT(item1.group.toLowerCase(), item2.group.toLowerCase());
                        if (g !== 0) {
                            return g;
                        }
                        return NATSORT(item1.text.toLowerCase(), item2.text.toLowerCase());
                    });
                    action.payload.features = features;
                }
            }
            return action.payload;
        case UPDATE_DATASET:
            if (action.payload.id === state.id) { // update name, description, url
                document.title = action.payload.name + ' - Cirro';
                return Object.assign(state, action.payload);
            }
            return state;
        default:
            return state;
    }
}

function datasetChoices(state = [], action) {
    switch (action.type) {
        case SET_DATASET_CHOICES:
            action.payload.sort(NATSORT);
            return action.payload;
        case SET_EMAIL:
            if (action.payload == null) {
                return [];
            }
            return state;
        case ADD_DATASET:
            state.push(action.payload);
            state.sort(NATSORT);
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
                    state.sort(NATSORT);
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

// set the selected embeddings, each embedding has name (str) e.g X_umap, dimensions (int), mode (str)
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


function selection(state = null, action) {
    switch (action.type) {
        case SET_SELECTION:
            return action.payload;
        case SET_DATASET:
            return null;
        default:
            return state;
    }
}

/**
 * Feature summary maps measure and dimension names to another object containing
 * categories and counts for dimensions, and statistics such as min, max for measures.
 * Features summaries are in the space of selected cells. Example:
 * featureSummary['louvain'] =
 *
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
            if (action.payload != null) {
                for (let key in action.payload) {
                    state[key] = action.payload[key];
                }
                return Object.assign({}, state);
            }
            return state;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

function serverInfo(state = {}, action) {
    switch (action.type) {
        case SET_SERVER_INFO:
            if (action.payload.ontology && action.payload.ontology.cellTypes) {
                action.payload.ontology.cellTypes.forEach(item => item.text = item.name);
                action.payload.ontology.cellTypes.sort((item1, item2) => NATSORT(item1.text, item2.text));
            }
            if (action.payload.datasetSelectorColumns == null) {
                action.payload.datasetSelectorColumns = [{field: 'name', label: 'Name', visible: true}, {
                    field: 'species',
                    label: 'Species',
                    visible: true
                }, {
                    field: 'title',
                    label: 'Title',
                    visible: true
                }];
            }
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


// one of embedding, distribution, composition, results
function tab(state = 'embedding', action) {
    switch (action.type) {
        case SET_TAB:
            return action.payload;
        case SET_DATASET:
            return 'embedding';
        case SET_DISTRIBUTION_DATA:
            return (state === 'distribution' && action.payload.length === 0) ? 'embedding' : state;
        case SET_SEARCH_TOKENS:
            return (state === 'composition' && action.payload.filter(item => item.type === FEATURE_TYPE.OBS_CAT).length < 2) ? 'embedding' : state;
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


function distributionPlotInterpolator(state = DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR_OBJ, action) {
    switch (action.type) {
        case SET_DISTRIBUTION_PLOT_INTERPOLATOR:
            return Object.assign({}, DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR_OBJ, action.payload);
        case RESTORE_VIEW:
            if (action.payload.distributionPlotInterpolator != null) {
                const interpolator = Object.assign({}, DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR_OBJ, action.payload.distributionPlotInterpolator);
                interpolator.value = getInterpolator(interpolator.name);
                return interpolator;
            }
            return state;
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


function jobResultId(state = null, action) {
    switch (action.type) {
        case SET_DATASET:
            return null;
        case SET_JOB_RESULT:
            return action.payload;
        default:
            return state;
    }
}


function jobResults(state = [], action) {
    switch (action.type) {
        case SET_DATASET:
            return [];
        case SET_JOB_RESULTS:
            return action.payload;
        default:
            return state;
    }
}

function distributionData(state = {}, action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            return Object.assign({}, state);
        case SET_DISTRIBUTION_DATA:
            return action.payload;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}

function selectedDistributionData(state = {}, action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            return Object.assign({}, state);
        case SET_SELECTED_DISTRIBUTION_DATA:
            return action.payload;
        case SET_DATASET:
            return {};
        default:
            return state;
    }
}


function datasetFilters(state = [], action) {
    switch (action.type) {
        case SET_DATASET:
            return [];
        case SET_DATASET_FILTERS:
            return action.payload;
        default:
            return state;
    }
}

function datasetViews(state = [], action) {
    switch (action.type) {
        case SET_DATASET:
            return [];
        case SET_DATASET_VIEWS:
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


// category -> originalValue -> {newValue, positiveMarkers, negativeMarkers, color}
function categoricalNames(state = {}, action) {
    switch (action.type) {
        case SET_CATEGORICAL_NAME:
            return action.payload;
        case UPDATE_CATEGORICAL_NAME:
            let category = state[action.payload.name];
            if (category === undefined) {
                category = {};
                state[action.payload.name] = category;
            }
            category[action.payload.originalValue] = action.payload;
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
            state.forEach((trace, index) => {
                if (trace.type === TRACE_TYPE_META_IMAGE) {
                    trace.categoryToStats = action.payload.size === 0 ? trace.fullCategoryToStats : createCategoryToStats(trace, action.payload);
                    updateTraceColors(trace);
                    state[index] = Object.assign({}, trace);
                }
            });
            return state.slice();
        case SET_DOMAIN:
            state.forEach((trace, index) => {
                if (trace.continuous && trace.name === action.payload.name) {
                    const summary = action.payload.summary;
                    const domain = trace.type === TRACE_TYPE_META_IMAGE ? [-3, 3] : [summary.min, summary.max];

                    if (trace.type === TRACE_TYPE_META_IMAGE) {
                        if (summary.customZMin != null && !isNaN(summary.customZMin)) {
                            domain[0] = summary.customZMin;
                        }
                        if (summary.customZMax != null && !isNaN(summary.customZMax)) {
                            domain[1] = summary.customZMax;
                        }
                    } else {
                        if (summary.customMin != null && !isNaN(summary.customMin)) {
                            domain[0] = summary.customMin;
                        }
                        if (summary.customMax != null && !isNaN(summary.customMax)) {
                            domain[1] = summary.customMax;
                        }
                    }

                    trace.colorScale.domain(domain);
                    updateTraceColors(trace);
                    state[index] = Object.assign({}, trace);
                }
            });
            return state.slice();
        case UPDATE_CATEGORICAL_COLOR:
            state.forEach((trace, index) => {
                if (!trace.continuous && trace.name === action.payload.name) {
                    const range = trace.colorScale.range();
                    range[trace.colorScale.domain().indexOf(action.payload.originalValue)] = action.payload.color;
                    trace.colorScale.range(range);
                    updateTraceColors(trace);
                    state[index] = Object.assign({}, trace);
                }
            });
            return state.slice();
        case SET_INTERPOLATOR:
            // update colors for existing continuous traces
            state.forEach((trace, index) => {
                if (trace.continuous && trace.featureType === action.payload.featureType) {
                    let domain = trace.colorScale.domain();
                    trace.colorScale = createColorScale(action.payload.value).domain(domain);
                    updateTraceColors(trace);
                    state[index] = Object.assign({}, trace);
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

export function unselectedPointSize(state = DEFAULT_POINT_SIZE, action) {
    switch (action.type) {
        case SET_UNSELECTED_POINT_SIZE:
            return action.payload;
        case RESTORE_VIEW:
            return action.payload.unselectedPointSize != null ? action.payload.unselectedPointSize : state;
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


// object with name, type, embeddingKey
export function activeFeature(state = {}, action) {
    switch (action.type) {
        case SET_DATASET:
            return null;
        case SET_ACTIVE_FEATURE:
            return action.payload;
        default:
            return state;
    }
}

// feature name (e.g. leiden) -> scroll position for restoring state
export function legendScrollPosition(state = {}, action) {
    switch (action.type) {
        case SET_DATASET:
            return {};
        case SET_LEGEND_SCROLL_POSITION:
            state[action.payload.name] = action.payload.value;
            return state;
        default:
            return state;
    }
}

// array of {label}
function tasks(state = [], action) {
    switch (action.type) {
        case ADD_TASK:
            return [...state, action.payload];
        case REMOVE_TASK:
            state.splice(state.indexOf(action.payload), 1);
            return state.slice();
        default:
            return state;
    }
}


function interpolator(state = DEFAULT_INTERPOLATORS, action) {
    switch (action.type) {
        case SET_INTERPOLATOR:
            const newValue = {};
            newValue[action.payload.featureType] = action.payload.value;
            return Object.assign({}, state, newValue);
        case RESTORE_VIEW:
            return action.payload.interpolator != null ? Object.assign({}, state, action.payload.interpolator) : state;
        default:
            return state;
    }
}

export default combineReducers({
    activeFeature,
    cachedData,
    categoricalNames,
    chartOptions,
    chartSize,
    combineDatasetFilters,
    dataset,
    datasetChoices,
    datasetFilter,
    datasetFilters,
    datasetViews,
    dialog,
    distributionData,
    distributionPlotOptions,
    distributionPlotInterpolator,
    email,
    embeddingData,
    embeddingLabels,
    embeddings,
    featureSummary,
    globalFeatureSummary,
    interpolator,
    jobResultId,
    jobResults,
    legendScrollPosition,
    tasks,
    loadingApp,
    markerOpacity,
    markers,
    message,
    panel,
    pointSize,
    savedDatasetState,
    searchTokens,
    selectedDistributionData,
    selection,
    serverInfo,
    tab,
    unselectedMarkerOpacity,
    unselectedPointSize,
    user
});

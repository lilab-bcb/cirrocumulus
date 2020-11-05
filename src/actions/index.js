import {scaleOrdinal, scaleSequential} from 'd3-scale';
import {schemeCategory10, schemePaired} from 'd3-scale-chromatic';
import {saveAs} from 'file-saver';
import {isEqual, uniqBy} from 'lodash';
import OpenSeadragon from 'openseadragon';
import CustomError from '../CustomError';
import {JsonDataset} from '../JsonDataset';
import {RestDataset} from '../RestDataset';
import {RestServerApi} from '../RestServerApi';


import {getPositions} from '../ThreeUtil';

import {
    addFeatureSetsToX,
    CATEGORY_20B,
    CATEGORY_20C,
    convertBinsToPoints,
    getFeatureSets,
    getInterpolator,
    indexSort,
    randomSeq,
    splitSearchTokens,
    updateTraceColors
} from '../util';

// export const API = 'http://localhost:5000/api';
export const API = 'api';

const authScopes = [
    'email',
    // 'profile',
    // 'https://www.googleapis.com/auth/userinfo.profile',
    // 'https://www.googleapis.com/auth/contacts.readonly',
    // 'https://www.googleapis.com/auth/devstorage.full_control',
];
export const SET_DOT_PLOT_INTERPOLATOR = 'SET_DOT_PLOT_INTERPOLATOR';
export const SET_CHART_OPTIONS = 'SET_CHART_OPTIONS';
export const SET_COMBINE_DATASET_FILTERS = 'SET_COMBINE_DATASET_FILTERS';
export const SET_DATASET_FILTERS = 'SET_DATASET_FILTERS'; // saved dataset filters
export const SET_UNSELECTED_MARKER_SIZE = 'SET_UNSELECTED_MARKER_SIZE';
export const SET_UNSELECTED_MARKER_SIZE_UI = 'SET_UNSELECTED_MARKER_SIZE_UI';
export const SET_PRIMARY_TRACE_KEY = 'SET_PRIMARY_TRACE_KEY';
export const SET_CHART_SIZE = 'SET_CHART_SIZE';
export const SET_PRIMARY_CHART_SIZE = 'SET_PRIMARY_CHART_SIZE';
export const SET_SERVER_INFO = "SET_SERVER_INFO";
export const SET_DATASET_FILTER = 'SET_DATASET_FILTER';
export const ADD_DATASET = 'ADD_DATASET';
export const DELETE_DATASET = 'DELETE_DATASET';
export const UPDATE_DATASET = 'UPDATE_DATASET';
export const SET_GLOBAL_FEATURE_SUMMARY = 'SET_GLOBAL_FEATURE_SUMMARY';
export const SET_SAVED_DATASET_STATE = 'SET_SAVED_DATASET_STATE';

export const SET_DOMAIN = 'SET_DOMAIN';
export const SET_CATEGORICAL_COLOR = 'SET_CATEGORICAL_COLOR';
export const SET_CATEGORICAL_NAME = 'SET_CATEGORICAL_NAME';
export const SET_MARKER_SIZE = 'SET_MARKER_SIZE';
export const SET_MARKER_OPACITY = 'SET_MARKER_OPACITY';

export const SET_EMBEDDING_CHART_SIZE = "SET_EMBEDDING_CHART_SIZE";

export const SET_UNSELECTED_MARKER_OPACITY = 'SET_UNSELECTED_MARKER_OPACITY';

// update chart

export const SET_SELECTION = 'SET_SELECTION';
export const SET_FEATURE_SUMMARY = 'SET_FEATURE_SUMMARY';
export const SET_SAVED_DATASET_FILTER = 'SET_SAVED_DATASET_FILTER';
export const SET_SEARCH_TOKENS = 'SET_SEARCH_TOKENS';

export const SET_SELECTED_EMBEDDING = 'SET_SELECTED_EMBEDDING';
export const SET_NUMBER_OF_BINS = 'SET_NUMBER_OF_BINS';
export const SET_BIN_VALUES = 'SET_BIN_VALUES';
export const SET_BIN_SUMMARY = 'SET_BIN_SUMMARY';
export const SET_MESSAGE = 'SET_MESSAGE';
export const SET_INTERPOLATOR = 'SET_INTERPOLATOR';
export const SET_POINT_SIZE = 'SET_POINT_SIZE';

// update ui
export const SET_EMAIL = 'SET_EMAIL';
export const SET_USER = 'SET_USER';
export const SET_DATASET = 'SET_DATASET';
export const SET_MARKERS = 'SET_MARKERS';
export const SET_DIALOG = 'SET_DIALOG';

export const MORE_OPTIONS_DIALOG = 'MORE_OPTIONS_DIALOG';
export const EDIT_DATASET_DIALOG = 'EDIT_DATASET_DIALOG';
export const IMPORT_DATASET_DIALOG = 'IMPORT_DATASET_DIALOG';
export const SAVE_DATASET_FILTER_DIALOG = 'SAVE_DATASET_FILTER_DIALOG';
export const SAVE_FEATURE_SET_DIALOG = 'SAVE_FEATURE_SET_DIALOG';
export const HELP_DIALOG = 'HELP_DIALOG';
export const DELETE_DATASET_DIALOG = 'DELETE_DATASET_DIALOG';

export const SET_DATASET_CHOICES = 'SET_DATASET_CHOICES';
export const RESTORE_VIEW = 'RESTORE_VIEW';

export const SET_DOT_PLOT_DATA = 'SET_DOT_PLOT_DATA';
export const SET_SELECTED_DOT_PLOT_DATA = 'SET_SELECTED_DOT_PLOT_DATA';

export const SET_DOT_PLOT_SORT_ORDER = 'SET_DOT_PLOT_SORT_ORDER';
export const SET_EMBEDDING_DATA = 'SET_EMBEDDING_DATA';

export const SET_LOADING = 'SET_LOADING';
export const SET_TAB = 'SET_TAB';

export const SET_LOADING_APP = 'LOADING_APP';

export const SET_NUMBER_OF_BINS_UI = 'SET_NUMBER_OF_BINS_UI';
export const SET_MARKER_SIZE_UI = 'SET_MARKER_SIZE_UI';
export const SET_MARKER_OPACITY_UI = 'SET_MARKER_OPACITY_UI';
export const SET_UNSELECTED_MARKER_OPACITY_UI = 'SET_UNSELECTED_MARKER_OPACITY_UI';

export function getEmbeddingKey(embedding) {
    let fullName = embedding.name + '_' + embedding.dimensions;
    if (embedding.bin || embedding.precomputed) {
        fullName = fullName + '_' + embedding.nbins + '_' + embedding.agg;
    }
    return fullName;
}

export function getTraceKey(traceInfo) {
    return traceInfo.name + '_' + getEmbeddingKey(traceInfo.embedding);
}

function getEmbeddingJson(embedding) {
    let value = {basis: embedding.name, ndim: embedding.dimensions};
    if (embedding.precomputed) {
        value.precomputed = true;
    }
    if (embedding.bin) {
        value.nbins = embedding.nbins;
        value.agg = embedding.agg;
    }

    return value;
}


function getUser() {
    return function (dispatch, getState) {
        getState().serverInfo.api.getUserPromise().then(user => dispatch(setUser(user)));
    };
}

export function initGapi() {
    return function (dispatch, getState) {
        dispatch(_setLoadingApp({loading: true, progress: 0}));
        const startTime = new Date().getTime();
        const approximateColdBootTime = 32;

        function loadingAppProgress() {
            if (getState().loadingApp.loading) {
                let elapsed = (new Date().getTime() - startTime) / 1000;
                let p = Math.min(100, 100 * (elapsed / approximateColdBootTime));
                dispatch(_setLoadingApp({loading: true, progress: p}));
                if (p < 100) {
                    window.setTimeout(loadingAppProgress, 300);
                }
            }
        }

        window.setTimeout(loadingAppProgress, 500);

        return fetch(API + '/server').then(result => result.json()).then(serverInfo => {
            serverInfo.api = new RestServerApi();
            if (serverInfo.clientId == null) {
                serverInfo.clientId = '';
            }
            if (serverInfo.fancy == null) {
                serverInfo.fancy = true;
            }
            dispatch(setServerInfo(serverInfo));
            if (serverInfo.clientId === '') { // no login required
                dispatch(_setLoadingApp({loading: false}));
                dispatch(_setEmail(serverInfo.mode === 'server' ? '' : null));
                if (serverInfo.mode === 'server') {
                    dispatch(setUser({importer: true}));
                }
                dispatch(listDatasets()).then(() => {
                    dispatch(_loadSavedView());
                });
            } else {
                console.log((new Date().getTime() - startTime) / 1000 + " startup time");
                let script = document.createElement('script');
                script.type = 'text/javascript';
                script.src = 'https://apis.google.com/js/api.js';
                script.onload = (e) => {
                    window.gapi.load('client:auth2', () => {
                        window.gapi.client.init({
                            clientId: serverInfo.clientId,
                            scope: authScopes.join(' '),
                        }).then(() => {
                            dispatch(_setLoadingApp({loading: false, progress: 0}));
                            dispatch(initLogin(true));
                        });
                    });
                };
                document.getElementsByTagName('head')[0].appendChild(script);
            }
        }).catch(err => {
            console.log(err);
            dispatch(_setLoadingApp({loading: false, progress: 0}));
            handleError(dispatch, 'Unable to load app. Please try again.');
        });
    };
}

export function logout() {
    return function (dispatch, getState) {
        window.gapi.auth2.getAuthInstance().signOut().then(() => {
            dispatch({type: SET_EMAIL, payload: null});
            dispatch(_setDataset(null));
        });
    };
}

export function login() {
    return function (dispatch, getState) {
        window.gapi.auth2.getAuthInstance().signIn().then(e => {
            dispatch(initLogin());
        });
    };
}

export function openDatasetFilter(filterId) {
    return function (dispatch, getState) {

        dispatch(_setLoading(true));

        getState().serverInfo.api.getDatasetFilterPromise(filterId, getState().dataset.id)
            .then(result => {
                result.id = filterId;
                let filterValue = JSON.parse(result.value);
                let datasetFilter = {};

                filterValue.filters.forEach(item => {
                    let filterOperation = item[1];
                    let field = item[0];
                    let slashIndex = field.indexOf('/');
                    if (slashIndex !== -1) {
                        field = field.substring(slashIndex + 1);
                    }
                    if (filterOperation === 'in') {
                        datasetFilter[field] = item[2];
                    } else {
                        datasetFilter[field] = {operation: item[1], value: item[2]};
                    }
                });
                dispatch(setDatasetFilter(datasetFilter));
                dispatch(handleFilterUpdated());


            }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve filter. Please try again.');
        });
    };
}

export function deleteDatasetFilter(filterId) {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        getState().serverInfo.api.deleteDatasetFilterPromise(filterId, getState().dataset.id)
            .then(() => {
                let datasetFilters = getState().datasetFilters;
                for (let i = 0; i < datasetFilters.length; i++) {
                    if (datasetFilters[i].id === filterId) {
                        datasetFilters.splice(i, 1);
                        break;
                    }
                }

                dispatch(setDatasetFilters(datasetFilters.slice()));
                dispatch(setMessage('Filter deleted'));

            }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete filter. Please try again.');
        });
    };
}

export function deleteFeatureSet(id) {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        getState().serverInfo.api.deleteFeatureSet(id, getState().dataset.id)
            .then(result => {

                let markers = getState().markers;
                for (let i = 0; i < markers.length; i++) {
                    if (markers[i].id === id) {
                        markers.splice(i, 1);
                        // remove from searchTokens
                        break;
                    }
                }
                dispatch(setMarkers(markers.slice()));
                dispatch(setMessage('Set deleted'));
            }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete set. Please try again.');
        });
    };
}

export function saveFeatureSet(payload) {
    return function (dispatch, getState) {
        const state = getState();
        const searchTokens = state.searchTokens;
        const splitTokens = splitSearchTokens(searchTokens);
        let features = splitTokens.X;
        const requestBody = {
            ds_id: state.dataset.id,
            name: payload.name,
            features: features,
            category: payload.category
        };

        dispatch(_setLoading(true));
        getState().serverInfo.api.upsertFeatureSet(requestBody, false)
            .then(result => {

                let markers = getState().markers;
                markers.push({
                    category: payload.category,
                    name: payload.name,
                    features: features,
                    id: result.id
                });
                dispatch(setMarkers(markers.slice()));
                dispatch(setMessage('Set saved'));
            }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to save set. Please try again.');
        });
    };
}

export function saveDatasetFilter(payload) {
    return function (dispatch, getState) {
        const state = getState();
        let filterValue = getFilterJson(state, true);

        let requestBody = {ds_id: state.dataset.id};
        if (payload.filterId != null) {
            requestBody.id = payload.filterId;
        }
        let savedDatasetFilter = {};

        let sendRequest = false;
        if (payload.name !== savedDatasetFilter.name) {
            requestBody.name = payload.name;
            sendRequest = true;
        }
        if (payload.notes !== savedDatasetFilter.notes) {
            requestBody.notes = payload.notes;
            sendRequest = true;
        }
        if (!isEqual(filterValue, savedDatasetFilter.value)) {
            requestBody.value = filterValue;
            sendRequest = true;
        }
        if (!sendRequest) {
            dispatch(setDialog(null));
            return;
        }
        for (let key in requestBody) {
            savedDatasetFilter[key] = requestBody[key];
        }
        dispatch(_setLoading(true));
        getState().serverInfo.api.upsertDatasetFilterPromise(requestBody, payload.filterId != null)
            .then(result => {
                requestBody.id = result.id;
                let datasetFilters = getState().datasetFilters;
                if (payload.filterId != null) {
                    for (let i = 0; i < datasetFilters.length; i++) {
                        if (datasetFilters[i].id === result.id) {
                            datasetFilters[i] = requestBody;
                            break;
                        }
                    }
                } else {
                    datasetFilters.push(requestBody);
                }
                dispatch(setDatasetFilters(datasetFilters.slice()));
                dispatch(setMessage('Filter saved'));
            }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to save filter. Please try again.');
        });
    };
}

export function removeDatasetFilter(filterKey) {
    return function (dispatch, getState) {
        if (filterKey == null) {
            // clear all
            dispatch(setDatasetFilter({}));
        } else {
            let datasetFilter = getState().datasetFilter;
            if (filterKey === 'selection') {
                for (let key in datasetFilter) {
                    if (datasetFilter[key].basis != null) {
                        delete datasetFilter[key];
                    }
                }
            } else {
                delete datasetFilter[filterKey];
            }
            dispatch(setDatasetFilter(Object.assign({}, datasetFilter)));
        }
        dispatch(handleFilterUpdated());
    };
}

export function getDatasetFilterArray(datasetFilter) {
    let filters = [];
    for (let key in datasetFilter) {
        // basis, path for brush filter
        const value = datasetFilter[key];
        let f = null;
        if (window.Array.isArray(value)) {
            f = [key, 'in', value];
        } else if (value.basis != null) {
            let points = value.points;
            if (points.length > 1) { // take "or"
                let allPoints = new Set();
                points.forEach(p => p.forEach(i => allPoints.add(i)));
                points = Array.from(allPoints);
            } else {
                points = points[0];
            }

            f = [getEmbeddingJson(value.basis), 'in', {points: points}];
        } else {
            if (value.operation !== '' && !isNaN(value.value) && value.value != null) {
                f = [key, value.operation, value.value];
            }
        }

        if (f != null) {
            filters.push(f);
        }
    }
    return filters;
}


function getFilterJson(state) {
    let filters = getDatasetFilterArray(state.datasetFilter);
    if (filters.length > 0) {
        const obs = state.dataset.obs;
        const obsCat = state.dataset.obsCat;
        for (let i = 0; i < filters.length; i++) {
            // add obs/ prefix
            if (obsCat.indexOf(filters[i][0]) !== -1 || obs.indexOf(filters[i][0]) !== -1) {
                filters[i][0] = 'obs/' + filters[i][0];
            }
        }
        return {filters: filters, combine: state.combineDatasetFilters};
    }
}


export function downloadSelectedIds() {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        let filter = getFilterJson(getState(), true);
        getState().dataset.api.getSelectedIdsPromise({
            filter: filter
        }).then(result => {
            const blob = new Blob([result.ids.join('\n')], {type: "text/plain;charset=utf-8"});
            saveAs(blob, "selection.txt");
        }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err);
        });
    };
}

export function exportDatasetFilters() {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        getState().serverInfo.api.exportDatasetFiltersPromise(getState().dataset.id).then(result => {
            if (result == null) {
                handleError(dispatch, 'Unable to export filters');
                return;
            }
            const blob = new Blob([result], {type: "text/plain;charset=utf-8"});
            saveAs(blob, "filters.csv");
        }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err);
        });
    };
}


export function setDotPlotSortOrder(payload) {
    return {type: SET_DOT_PLOT_SORT_ORDER, payload: payload};
}

export function setChartOptions(payload) {
    return {type: SET_CHART_OPTIONS, payload: payload};
}


function _setCombineDatasetFilters(payload) {
    return {type: SET_COMBINE_DATASET_FILTERS, payload: payload};
}

export function setCombineDatasetFilters(payload) {
    return function (dispatch, getState) {
        dispatch(_setCombineDatasetFilters(payload));
        dispatch(handleFilterUpdated());
    };
}

export function setChartSize(payload) {
    return {type: SET_CHART_SIZE, payload: payload};
}


export function setPrimaryChartSize(payload) {
    return {type: SET_PRIMARY_CHART_SIZE, payload: payload};
}

export function setPrimaryTraceKey(payload) {
    return {type: SET_PRIMARY_TRACE_KEY, payload: payload};
}

function setGlobalFeatureSummary(payload) {
    return {type: SET_GLOBAL_FEATURE_SUMMARY, payload: payload};
}

function setDatasetFilter(payload) {
    return {type: SET_DATASET_FILTER, payload: payload};
}

function setSelection(payload) {
    return {type: SET_SELECTION, payload: payload};
}

function setFeatureSummary(payload) {
    return {type: SET_FEATURE_SUMMARY, payload: payload};
}


function handleFilterUpdated() {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        // whenever filter is updated, we need to get selection statistics
        const state = getState();

        const searchTokens = splitSearchTokens(state.searchTokens);
        let filter = getFilterJson(state);
        addFeatureSetsToX(getFeatureSets(state.markers, searchTokens.featureSets), searchTokens.X);
        let q = {
            selection: {
                measures: searchTokens.X.concat(searchTokens.obs.map(item => 'obs/' + item)),
                dimensions: searchTokens.obsCat
            }
        };

        if (filter) {
            q.selection.filter = filter;
        }
        const isCurrentSelectionEmpty = state.selection.chart == null || Object.keys(state.selection.chart).length === 0;

        if (filter == null) {
            state.dotPlotData.forEach(data => {
                data.selection = null;
            });
            dispatch(_setDotPlotData(state.dotPlotData.slice()));
            if (!isCurrentSelectionEmpty) {
                dispatch(setSelection({chart: {}}));
            }
            dispatch(_setSelectedDotPlotData([]));
            dispatch(setFeatureSummary({}));
            dispatch(_setLoading(false));
            return;
        }
        let selectedEmbeddings = state.embeddings;
        q.selection.embeddings = selectedEmbeddings.map(e => {
            return getEmbeddingJson(e);
        });

        getState().dataset.api.getDataPromise(q).then(result => {
            dispatch(handleSelectionResult(result.selection, true));
        }).catch(err => {
            handleError(dispatch, err);
        }).finally(() => dispatch(_setLoading(false)));
    };
}


export function handleBrushFilterUpdated(payload) {
    return function (dispatch, getState) {
        const name = payload.name; // full basis name
        const value = payload.value;  // value has basis and points
        const clear = payload.clear;
        let datasetFilter = getState().datasetFilter;

        let update = true;
        if (value == null) { // remove
            update = datasetFilter[name] != null;
            delete datasetFilter[name];
        } else {
            if (clear) {
                datasetFilter[name] = {basis: value.basis, points: [value.points]};
            } else {
                const prior = datasetFilter[name];
                if (prior != null) {
                    datasetFilter[name].points.push(value.points);
                } else {
                    datasetFilter[name] = {basis: value.basis, points: [value.points]};
                }
            }
        }
        if (update) {
            dispatch(setDatasetFilter(Object.assign({}, datasetFilter)));
            dispatch(handleFilterUpdated());
        }
    };
}


export function handleMeasureFilterUpdated(payload) {
    return function (dispatch, getState) {
        const name = payload.name;
        const operation = payload.operation;
        const value = payload.value;
        let update = payload.update;

        let datasetFilter = getState().datasetFilter;
        let filter = datasetFilter[name];

        if (filter == null) {
            filter = {operation: '>', value: NaN};
            datasetFilter[name] = filter;
        }
        if (update) {
            if (value != null) {
                update = (!isNaN(value) ? value !== filter.value : isNaN(value) !== isNaN(filter.value));
            }
            if (operation != null) {
                update = update || (operation !== filter.operation && !isNaN(filter.value));
            }
        }
        if (operation != null) {
            filter.operation = operation;
        }
        if (value != null) {
            filter.value = value;
        }

        dispatch(setDatasetFilter(Object.assign({}, datasetFilter)));
        if (update) {
            dispatch(handleFilterUpdated());
        }
    };
}

export function handleColorChange(payload) {
    return {type: SET_CATEGORICAL_COLOR, payload: payload};
}

export function handleDomainChange(payload) {
    return {type: SET_DOMAIN, payload: payload};
}

function _handleCategoricalNameChange(payload) {
    return {type: SET_CATEGORICAL_NAME, payload: payload};

}

export function handleCategoricalNameChange(payload) {
    return function (dispatch, getState) {
        const dataset_id = getState().dataset.id;
        const p = getState().serverInfo.fancy ? getState().serverInfo.api.setCategoryNamePromise({
            c: payload.name,
            o: payload.oldValue,
            n: payload.value,
            id: dataset_id
        }) : Promise.resolve();
        p.then(() => dispatch(_handleCategoricalNameChange(payload)));
    };
}

export function handleDimensionFilterUpdated(payload) {
    return function (dispatch, getState) {
        let name = payload.name;
        let value = payload.value;
        let shiftKey = payload.shiftKey;
        let metaKey = payload.metaKey;
        let datasetFilter = getState().datasetFilter;
        let embeddingData = getState().embeddingData;
        let categories;
        for (let i = 0; i < embeddingData.length; i++) {
            if (embeddingData[i].name === name) {
                categories = embeddingData[i].colorScale.domain();
                break;
            }
        }
        let selectedValues = datasetFilter[name];
        if (selectedValues == null) {
            selectedValues = [];
            datasetFilter[name] = selectedValues;
        }

        if (shiftKey && selectedValues.length > 0) {
            // add last click to current
            let lastIndex = categories.indexOf(selectedValues[selectedValues.length - 1]);
            let currentIndex = categories.indexOf(value);
            // put clicked category at end of array
            if (currentIndex > lastIndex) {
                for (let i = lastIndex; i <= currentIndex; i++) {
                    let index = selectedValues.indexOf(value);
                    if (index !== -1) {
                        selectedValues.splice(index, 1);
                    }
                    selectedValues.push(categories[i]);
                }
            } else {
                for (let i = lastIndex; i >= currentIndex; i--) {
                    let index = selectedValues.indexOf(value);
                    if (index !== -1) {
                        selectedValues.splice(index, 1);
                    }
                    selectedValues.push(categories[i]);
                }
            }
        } else {
            let selectedIndex = selectedValues.indexOf(value);
            if (!metaKey) { // clear and toggle current
                selectedValues = [];
                datasetFilter[name] = selectedValues;
            }
            if (selectedIndex !== -1) { // exists, remove
                selectedValues.splice(selectedIndex, 1);
                if (selectedValues.length === 0) {
                    delete datasetFilter[name];
                }
            } else {
                selectedValues.push(value);
            }
        }
        dispatch(setDatasetFilter(datasetFilter));
        dispatch(handleFilterUpdated());
    };
}


export function restoreView(payload) {
    return {type: RESTORE_VIEW, payload: payload};
}

export function setNumberOfBinsUI(payload) {
    return {type: SET_NUMBER_OF_BINS_UI, payload: payload};
}

export function setPointSize(payload) {
    return {type: SET_POINT_SIZE, payload: payload};
}

export function setMarkerSize(payload) {
    return {type: SET_MARKER_SIZE, payload: payload};
}


export function setUser(payload) {
    return {type: SET_USER, payload: payload};
}

export function setMarkerOpacityUI(payload) {
    return {type: SET_MARKER_OPACITY_UI, payload: payload};
}

export function setUnselectedMarkerOpacityUI(payload) {
    return {type: SET_UNSELECTED_MARKER_OPACITY_UI, payload: payload};
}

export function setServerInfo(payload) {
    return {type: SET_SERVER_INFO, payload: payload};
}

function loadDefaultDatasetEmbedding() {
    return function (dispatch, getState) {
        const dataset = getState().dataset;
        const embeddingNames = dataset.embeddings.map(e => e.name);
        let priority = {'tissue_hires': 1, 'X_fle': 2, 'X_umap': 3, 'X_tsne': 4, 'X_fitsne': 5};

        embeddingNames.sort((a, b) => {
            a = priority[a] || Number.MAX_VALUE;
            b = priority[b] || Number.MAX_VALUE;
            return a - b;
        });

        if (embeddingNames.length > 0) {
            dispatch(setSelectedEmbedding([dataset.embeddings[dataset.embeddings.map(e => e.name).indexOf(embeddingNames[0])]]));
        }
    };
}

function loadDefaultDataset() {
    return function (dispatch, getState) {
        if (getState().dataset == null && getState().datasetChoices.length === 1) {
            dispatch(setDataset(getState().datasetChoices[0].id));
        }
    };
}

function restoreSavedView(savedView) {
    return function (dispatch, getState) {
        if (savedView.colorScheme != null) {
            let interpolator = getInterpolator(savedView.colorScheme);
            if (interpolator != null) {
                savedView.colorScheme = {
                    name: savedView.colorScheme,
                    value: interpolator,
                };
            }
        }

        if (savedView.datasetFilter != null) {
            for (let key in savedView.datasetFilter) {
                let value = savedView.datasetFilter[key];
                if (!window.Array.isArray(value)) {
                    value.uiValue = value.value;
                }
            }
        } else {
            savedView.datasetFilter = {};
        }

        dispatch(_setLoading(true));
        let datasetPromise;
        if (savedView.dataset != null) {
            datasetPromise = dispatch(setDataset(savedView.dataset, false, false));
        } else {
            datasetPromise = Promise.resolve();
        }
        datasetPromise
            .then(() => {
                let dataset = getState().dataset;
                if (savedView.embeddings) {
                    let names = dataset.embeddings.map(e => getEmbeddingKey(e));
                    let embeddings = [];
                    savedView.embeddings.forEach(embedding => {
                        let index = names.indexOf(getEmbeddingKey(embedding));
                        if (index !== -1) {
                            embeddings.push(dataset.embeddings[index]);
                        }
                    });
                    savedView.embeddings = embeddings;
                    if (savedView.camera != null) {
                        savedView.chartOptions.camera = savedView.camera;
                    }
                }
            })
            .then(() => dispatch(setDatasetFilter(savedView.datasetFilter)))
            .then(() => dispatch(restoreView(savedView)))
            .then(() => dispatch(_updateCharts()))
            .then(() => dispatch(handleFilterUpdated()))
            .then(() => {
                if (savedView.sort != null) {
                    for (let name in savedView.sort) {
                        dispatch(setDotPlotSortOrder({name: name, value: savedView.sort[name]}));
                    }
                }

            })
            .finally(() => dispatch(_setLoading(false)))
            .catch(err => {
                console.log(err);
                dispatch(setMessage('Unable to restore saved view.'));
                dispatch(loadDefaultDataset());
            });
    };
}

function _loadSavedView() {
    return function (dispatch, getState) {
        let savedView = {dataset: null};
        let q = window.location.search.substring(3);
        if (q.length > 0) {
            try {
                savedView = JSON.parse(window.decodeURIComponent(q));
            } catch (err) {
                return dispatch(setMessage('Unable to restore saved view.'));
            }
        }

        if (savedView.dataset != null) {
            dispatch(restoreSavedView(savedView));
        } else {
            dispatch(loadDefaultDataset());
        }
    };
}

export function initLogin(loadSavedView) {
    return function (dispatch, getState) {
        let isSignedIn = window.gapi.auth2.getAuthInstance().isSignedIn.get();
        if (isSignedIn) {
            let googleUser = window.gapi.auth2.getAuthInstance().currentUser.get();
            let email = googleUser.getBasicProfile().getEmail();
            dispatch(_setEmail(email));
            dispatch(getUser());
            dispatch(listDatasets()).then(() => {
                if (loadSavedView) {
                    dispatch(_loadSavedView());
                }
            });

        }
    };
}

export function getAccessToken() {
    return window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().access_token;
}

export function getIdToken() {
    return typeof window.gapi !== 'undefined' ? window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().id_token : '';
}

export function saveDataset(payload) {
    return function (dispatch, getState) {

        const name = payload.name;
        const url = payload.url;

        if (name === '' || url === '') {
            return;
        }
        const readers = payload.readers;
        const description = payload.description;
        // let bucket = url.substring('gs://'.length);
        // let slash = bucket.indexOf('/');
        // let object = encodeURIComponent(bucket.substring(slash + 1));
        // bucket = encodeURIComponent(bucket.substring(0, slash));

        dispatch(_setLoading(true));
        const isEdit = payload.dataset != null;
        const serverInfo = getState().serverInfo;
        const updatePermissions = !isEdit || url !== payload.dataset.url;
        // if (updatePermissions) {

        // updateDatasetPermissionPromise = fetch(
        //     'https://www.googleapis.com/storage/v1/b/' + bucket + '/o/' + object,
        //     {
        //         body: JSON.stringify({
        //             'acl': [{
        //                 'entity': 'user-' + serverEmail,
        //                 'role': 'READER'
        //             }]
        //         }),
        //         method: 'PATCH',
        //         headers: {'Authorization': 'Bearer ' + getAccessToken(), 'Content-Type': 'application/json'},
        //     });

        // } else {
        let updateDatasetPermissionPromise = Promise.resolve({ok: true});
        // }
        updateDatasetPermissionPromise.then(permissionsResponse => {

            // if (!permissionsResponse.ok) {
            //     dispatch(setMessage('Unable to set dataset read permissions. Please ensure that you are the dataset owner or manually add ' + serverEmail + ' as a reader.'));
            // }
        }).then(() => getState().serverInfo.api.upsertDatasetPromise(payload.dataset ? payload.dataset.id : null,
            {name: name, readers: readers, description: description, url: url})).then(importDatasetResult => {
            if (isEdit) {
                dispatch(updateDataset({
                    name: name,
                    id: importDatasetResult.id,
                    description: description,
                    owner: true
                }));
                dispatch(setMessage('Dataset updated'));
            } else {
                dispatch(_addDataset({name: name, id: importDatasetResult.id, description: description, owner: true}));
                if (serverInfo.email) {
                    dispatch(setMessage(updatePermissions ? 'Please ensure that ' + serverInfo.email + ' is a "Storage Object Viewer"' : 'Dataset created'));
                }
            }

        }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err);
        });
    };
}

export function deleteDataset(payload) {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        getState().serverInfo.api.deleteDatasetPromise(payload.dataset.id).then(() => {
            dispatch(_setDataset(null));
            dispatch(_deleteDataset({id: payload.dataset.id}));
            dispatch(setDialog(null));

            dispatch(setMessage('Dataset deleted'));
        }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete dataset. Please try again.');
        });
    };
}

function _addDataset(payload) {
    return {type: ADD_DATASET, payload: payload};
}

function _deleteDataset(payload) {
    return {type: DELETE_DATASET, payload: payload};
}

export function updateDataset(payload) {
    return {type: UPDATE_DATASET, payload: payload};
}

export function setMessage(payload) {
    return {type: SET_MESSAGE, payload: payload};
}

export function setInterpolator(payload) {
    return {type: SET_INTERPOLATOR, payload: payload};
}
export function setDotPlotInterpolator(payload) {
    return {type: SET_DOT_PLOT_INTERPOLATOR, payload: payload};
}

export function setEmbeddingData(payload) {
    return {type: SET_EMBEDDING_DATA, payload: payload};
}

export function setDialog(payload) {
    return {type: SET_DIALOG, payload: payload};
}

export function setMarkers(payload) {
    return {type: SET_MARKERS, payload: payload};
}

function _setDatasetChoices(payload) {
    return {type: SET_DATASET_CHOICES, payload: payload};
}


function _setLoading(payload) {
    return {type: SET_LOADING, payload: payload};
}

export function setTab(payload) {
    return {type: SET_TAB, payload: payload};
}

function _setLoadingApp(payload) {
    return {type: SET_LOADING_APP, payload: payload};
}

function _setDataset(payload) {
    return {type: SET_DATASET, payload: payload};
}


export function setMarkerOpacity(payload) {
    return {type: SET_MARKER_OPACITY, payload: payload};
}

export function setUnselectedMarkerOpacity(payload) {
    return {type: SET_UNSELECTED_MARKER_OPACITY, payload: payload};
}


function _setDotPlotData(payload) {
    return {type: SET_DOT_PLOT_DATA, payload: payload};
}

function _setSelectedDotPlotData(payload) {
    return {type: SET_SELECTED_DOT_PLOT_DATA, payload: payload};
}


function _setEmail(payload) {
    return {type: SET_EMAIL, payload: payload};
}


/**
 *
 * @param value
 * @param type one of X, obs, featureSet
 */
export function setSearchTokens(value, type) {
    return function (dispatch, getState) {
        const state = getState();
        let searchTokens = state.searchTokens;
        // keep all other types
        let removeType = [type];
        if (type === 'obs') {
            removeType.push('obsCat');
        }
        searchTokens = searchTokens.filter(item => removeType.indexOf(item.type) === -1);
        if (type === 'X' || type === 'featureSet') {
            searchTokens = searchTokens.concat(value.map(item => {
                return {value: item, type: type};
            }));
        } else if (type === 'obs') {
            const obsCat = state.dataset.obsCat;
            value.forEach(val => {
                const type = obsCat.indexOf(val) !== -1 ? 'obsCat' : 'obs';
                searchTokens.push({value: val, type: type});
            });
        }
        dispatch({type: SET_SEARCH_TOKENS, payload: searchTokens.slice()});
        dispatch(_updateCharts());
    };
}


export function setSelectedEmbedding(payload) {
    return function (dispatch, getState) {
        let prior = getState().embeddings;
        dispatch({type: SET_SELECTED_EMBEDDING, payload: payload});
        dispatch(_updateCharts(err => {
            dispatch({type: SET_SELECTED_EMBEDDING, payload: prior});
        }));
    };
}


export function setNumberOfBins(payload) {
    return function (dispatch, getState) {
        if (getState().numberOfBins !== payload) {
            dispatch({type: SET_NUMBER_OF_BINS, payload: payload});
        }
    };
}

export function setBinSummary(payload) {
    return function (dispatch, getState) {
        dispatch({type: SET_BIN_SUMMARY, payload: payload});

    };
}

export function setBinValues(payload) {
    return {type: SET_BIN_VALUES, payload: payload};
}


export function setSavedDatasetState(payload) {
    return {type: SET_SAVED_DATASET_STATE, payload: payload};
}

function setDatasetFilters(payload) {
    return {type: SET_DATASET_FILTERS, payload: payload};
}


export function setDataset(id, loadDefaultView = true, setLoading = true) {
    return function (dispatch, getState) {
        if (setLoading) {
            dispatch(_setLoading(true));
        }
        let savedDatasetState = getState().savedDatasetState[id];
        const datasetChoices = getState().datasetChoices;
        let selectedChoice = null; // has id, owner, name
        for (let i = 0; i < datasetChoices.length; i++) {
            if (datasetChoices[i].id === id) {
                selectedChoice = datasetChoices[i];
                break;
            }
        }
        if (selectedChoice == null) {
            dispatch(_setLoading(false));
            dispatch(setMessage('Unable to find dataset'));
            return Promise.reject('Unable to find dataset');
        }
        // force re-render selected dataset dropdown
        let dataset = Object.assign({}, selectedChoice);
        dataset.id = id;
        dataset.embeddings = [];
        dataset.features = [];
        dataset.obs = [];
        dataset.obsCat = [];
        dispatch(_setDataset(dataset));
        let categoryNameResults;
        let datasetFilters = [];
        let newDataset;

        function onPromisesComplete() {
            newDataset = Object.assign(dataset, newDataset);
            if (newDataset.embeddings) {
                for (let i = 0; i < newDataset.embeddings.length; i++) {
                    if (newDataset.embeddings[i].nbins != null) {
                        newDataset.embeddings[i].bin = true;
                        newDataset.embeddings[i].precomputed = true;
                    } else { // set default values
                        newDataset.embeddings[i].nbins = 500;
                        newDataset.embeddings[i]._nbins = 500;
                        newDataset.embeddings[i].bin = false;
                        newDataset.embeddings[i].agg = 'max';
                    }
                }
            }
            newDataset.api = dataset.api;
            newDataset.features = newDataset.var;
            newDataset.features.sort((a, b) => {
                a = a.toLowerCase();
                b = b.toLowerCase();
                const aIsDigit = a[0] >= '0' && a[0] <= '9';
                const bIsDigit = b[0] >= '0' && b[0] <= '9';
                if (aIsDigit) {
                    a = 'zzzzzz' + a;
                }
                if (bIsDigit) {
                    b = 'zzzzzz' + b;
                }
                // put features that start with a number last
                return (a < b ? -1 : (a === b ? 0 : 1));
            });
            newDataset.id = id;
            dispatch(_setDataset(newDataset));

            if (categoryNameResults != null) {
                categoryNameResults.forEach(result => {
                    dispatch(_handleCategoricalNameChange({
                        name: result.category,
                        oldValue: result.original,
                        value: result.new
                    }));
                });
            }
            dispatch(setDatasetFilters(datasetFilters));
            if (savedDatasetState) {

                dispatch(restoreSavedView(savedDatasetState));
            } else if (loadDefaultView) {
                dispatch(loadDefaultDatasetEmbedding());
            }
        }

        let url = dataset.url || '';
        const isByteRangeDataset = url.startsWith("http") || url.startsWith('//');
        if (isByteRangeDataset) {
            dataset.api = new JsonDataset();
        } else {
            dataset.api = new RestDataset();
        }
        const initPromise = dataset.api.init(id, dataset.url);

        return initPromise.then(() => {
            let promises = [];
            if (!isByteRangeDataset) {
                const categoriesRenamePromise = getState().serverInfo.api.getCategoryNamesPromise(dataset.id).then(results => {
                    categoryNameResults = results;
                });

                const filtersPromise = getState().serverInfo.api.getFiltersPromise(dataset.id).then(results => {
                    datasetFilters = results;
                });
                promises.push(categoriesRenamePromise);
                promises.push(filtersPromise);
            }
            const schemaPromise = dataset.api.getSchemaPromise().then(result => {
                newDataset = result;
            });
            promises.push(schemaPromise);
            return promises;
        }).then(promises => Promise.all(promises)).then(() => onPromisesComplete()).finally(() => {
            if (setLoading) {
                dispatch(_setLoading(false));
            }
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve dataset. Please try again.');
        });

    };
}


function handleSelectionResult(selectionResult, clear) {
    return function (dispatch, getState) {
        const state = getState();
        if (selectionResult) {
            const coordinates = selectionResult.coordinates;
            if (coordinates != null) {
                let chartSelection = {};
                for (let key in coordinates) {
                    const selectedIndicesOrBins = coordinates[key].indices_or_bins;
                    const embedding = state.embeddings[state.embeddings.map(e => getEmbeddingKey(e)).indexOf(key)];
                    if (embedding == null) {
                        console.log(key + ' missing');
                        continue;
                    }
                    let selectedPoints = selectedIndicesOrBins;
                    if (embedding.bin) {
                        // find embedding data with matching layout
                        for (let i = 0; i < state.embeddingData.length; i++) {
                            const embedding = state.embeddingData[i].embedding;
                            let fullName = getEmbeddingKey(embedding);
                            if (fullName === key) {
                                const traceBins = state.embeddingData[i].bins;
                                if (traceBins != null) {
                                    selectedPoints = convertBinsToPoints(traceBins, selectedIndicesOrBins);
                                    break;
                                }
                            }
                        }

                    }

                    chartSelection[key] = {
                        userPoints: new Set(selectedPoints),
                        points: selectedIndicesOrBins
                    };
                }
                dispatch(setSelection({
                    count: selectionResult.count,
                    chart: chartSelection
                }));
            } else {
                const isCurrentSelectionEmpty = state.selection.chart == null || Object.keys(state.selection.chart).length === 0;
                if (clear && !isCurrentSelectionEmpty) {
                    dispatch(setSelection({count: selectionResult.count, chart: {}}));
                }
            }

            // userPoints are in chart space, points are in server space, count is total number of cells selected
            if (selectionResult.summary) {
                // merge or clear selection
                let selectionSummary = clear ? selectionResult.summary : Object.assign(getState().featureSummary, selectionResult.summary);
                dispatch(setFeatureSummary(Object.assign({}, selectionSummary)));
            }
            if (selectionResult.dotplot) {
                let selectedDotPlotData = state.selectedDotPlotData;
                if (clear) {
                    selectedDotPlotData = [];
                }
                selectedDotPlotData = updateDotPlotData(selectionResult.dotplot, selectedDotPlotData, splitSearchTokens(state.searchTokens));
                dispatch(_setSelectedDotPlotData(selectedDotPlotData));
            }

        }
    };

}

function updateDotPlotData(newDotplotData, dotPlotData, searchTokens) {
    let categoryKeys = [searchTokens.obsCat.join('-')];
    let dotPlotKeys = {}; // category-feature

    categoryKeys.forEach(category => {
        searchTokens.X.forEach(feature => {
            if (feature !== '__count') {
                let key = category + '-' + feature;
                dotPlotKeys[key] = true;

            }
        });
    });


    // merge with existing data
    if (newDotplotData) {
        newDotplotData.forEach(categoryItem => {
            let index = -1;
            for (let i = 0; i < dotPlotData.length; i++) {
                if (dotPlotData[i].name === categoryItem.name) {
                    index = i;
                    break;
                }
            }
            if (index === -1) {
                dotPlotData.push(categoryItem);
            } else {
                dotPlotData[index].values = categoryItem.values.concat(dotPlotData[index].values);
                dotPlotData[index].values = uniqBy(dotPlotData[index].values, (value => value.name));
                dotPlotData[index] = Object.assign({}, dotPlotData[index]);
            }
        });
    }

    // keep active categories only

    dotPlotData = dotPlotData.filter(categoryItem => categoryKeys.indexOf(categoryItem.name) !== -1);
    dotPlotData.forEach((categoryItem, categoryIndex) => {
        const nfeatures = categoryItem.values.length;
        categoryItem.values = categoryItem.values.filter(feature => dotPlotKeys[categoryItem.name + '-' + feature.name]);
        if (nfeatures !== categoryItem.values.length) {
            dotPlotData[categoryIndex] = Object.assign({}, categoryItem);
        }
    });
    dotPlotData = dotPlotData.filter(categoryItem => categoryItem.values.length > 0);
    // sort features
    let featureSortOrder = {};
    searchTokens.X.forEach((name, index) => {
        featureSortOrder[name] = index;
    });
    dotPlotData.forEach((categoryItem) => {
        categoryItem.values.sort((a, b) => {
            a = featureSortOrder[a.name];
            b = featureSortOrder[b.name];
            return a - b;
        });

    });
    return dotPlotData;
}

function _updateCharts(onError) {
    return function (dispatch, getState) {
        const state = getState();
        if (state.dataset == null) {
            return;
        }
        dispatch(_setLoading(true));

        const searchTokens = splitSearchTokens(state.searchTokens);
        if (searchTokens.featureSets.length > 0) {
            addFeatureSetsToX(getFeatureSets(state.markers, searchTokens.featureSets), searchTokens.X);
        }
        let dotplot = (searchTokens.X.length > 0 || searchTokens.featureSets.length > 0) && searchTokens.obsCat.length > 0;
        if (searchTokens.X.length === 0 && searchTokens.obsCat.length === 0 && searchTokens.obs.length === 0
            && searchTokens.featureSets.length === 0) {
            searchTokens.X = ['__count'];
        }
        const selectionStats = state.featureSummary;
        const embeddings = state.embeddings;
        let dotPlotData = state.dotPlotData;
        let selectedDotPlotData = state.selectedDotPlotData;
        let embeddingData = state.embeddingData;
        const globalFeatureSummary = state.globalFeatureSummary;
        let embeddingToVisibleFeatures = {};
        embeddings.forEach(selectedEmbedding => {
            let embeddingKey = getEmbeddingKey(selectedEmbedding);
            let features = {};
            searchTokens.X.concat(searchTokens.obs).concat(searchTokens.obsCat).forEach(feature => {
                features[feature] = true;
            });
            embeddingToVisibleFeatures[embeddingKey] = features;
        });

        // set active flag on cached data
        embeddingData.forEach(traceInfo => {
            const embeddingKey = getEmbeddingKey(traceInfo.embedding);
            let visibleFeatures = embeddingToVisibleFeatures[embeddingKey] || {};
            let active = visibleFeatures[traceInfo.name];
            if (active) {
                traceInfo.date = new Date();
            }
            // data is cached so delete from visibleFeatures
            delete visibleFeatures[traceInfo.name];
            traceInfo.active = active;
        });

        let dotplotMeasuresCacheMiss = new Set();
        let dotPlotCategories = searchTokens.obsCat.slice();
        let dotPlotCategoryKeys = [dotPlotCategories.join('-')];
        let dotPlotKeys = {}; // category-feature
        let selectedMeasuresCacheMiss = [];
        if (dotplot) {
            let cachedDotPlotKeys = {};
            dotPlotData.forEach(dotPlotDataItem => {
                dotPlotDataItem.values.forEach(datum => {
                    cachedDotPlotKeys[dotPlotDataItem.name + '-' + datum.name] = true;
                });
            });
            let cachedSelectedDotPlotMeasures = {};
            selectedDotPlotData.forEach(dotPlotDataItem => {
                dotPlotDataItem.values.forEach(datum => {
                    cachedSelectedDotPlotMeasures[datum.name] = true;
                });
            });

            searchTokens.X.forEach(feature => {
                if (feature !== '__count' && cachedSelectedDotPlotMeasures[feature] == null) {
                    selectedMeasuresCacheMiss.push(feature);
                }
            });
            searchTokens.obs.forEach(feature => {
                if (selectionStats[feature] == null) {
                    selectedMeasuresCacheMiss.push('obs/' + feature);
                }
            });
            dotPlotCategoryKeys.forEach(category => {
                let added = false;
                searchTokens.X.forEach(feature => {
                    if (feature !== '__count') {
                        let key = category + '-' + feature;
                        dotPlotKeys[key] = true;
                        if (cachedDotPlotKeys[key] == null) {
                            dotplotMeasuresCacheMiss.add(feature);
                            added = true;
                        }
                    }
                });
            });
        }

        let q = {};
        let embeddingsUsed = [];
        for (let embeddingName in embeddingToVisibleFeatures) {
            let visibleFeatures = embeddingToVisibleFeatures[embeddingName];
            let embedding = embeddings[embeddings.map(e => getEmbeddingKey(e)).indexOf(embeddingName)];
            let measures = [];
            let dimensions = [];
            for (let feature in visibleFeatures) {
                if (searchTokens.obs.indexOf(feature) !== -1) { // continuous obs
                    measures.push('obs/' + feature);
                } else if (searchTokens.obsCat.indexOf(feature) !== -1) { // categorical obs
                    dimensions.push(feature);
                } else {
                    measures.push(feature); // X
                }
            }
            if (measures.length > 0 || dimensions.length > 0) {
                embeddingsUsed.push(embedding);
                if (!q.embedding) {
                    q.embedding = [];
                }
                q.embedding.push(Object.assign({
                    measures: measures,
                    dimensions: dimensions
                }, getEmbeddingJson(embedding)));
            }
        }

        let globalFeatureSummaryXCacheMiss = [];
        let globalFeatureSummaryObsContinuousCacheMiss = [];
        let globalFeatureSummaryDimensionsCacheMiss = [];
        searchTokens.X.forEach(feature => {
            if (feature !== '__count') {
                if (globalFeatureSummary[feature] == null) {
                    globalFeatureSummaryXCacheMiss.push(feature);
                }
            }
        });
        searchTokens.obs.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryObsContinuousCacheMiss.push('obs/' + feature);
            }
        });
        searchTokens.obsCat.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryDimensionsCacheMiss.push(feature);
            }
        });


        if (globalFeatureSummaryXCacheMiss.length > 0 || globalFeatureSummaryObsContinuousCacheMiss.length > 0 || globalFeatureSummaryDimensionsCacheMiss.length > 0) {
            q.stats = {
                measures: globalFeatureSummaryXCacheMiss.concat(globalFeatureSummaryObsContinuousCacheMiss),
                dimensions: globalFeatureSummaryDimensionsCacheMiss
            };
        }

        if (dotPlotCategories.length > 0 && dotplotMeasuresCacheMiss.size > 0) {
            q.groupedStats = {
                measures: Array.from(dotplotMeasuresCacheMiss),
                dimensions: [dotPlotCategories]
            };
        }

        const filterJson = getFilterJson(state);
        // TODO update selection in new embedding space
        if (filterJson != null && (globalFeatureSummaryXCacheMiss.length > 0 || globalFeatureSummaryDimensionsCacheMiss.length > 0 || selectedMeasuresCacheMiss.length > 0)) {
            q.selection = {
                filter: filterJson,
                measures: globalFeatureSummaryDimensionsCacheMiss.length > 0 ? searchTokens.X : selectedMeasuresCacheMiss,
                dimensions: searchTokens.obsCat
            };
        }
        let dataPromise = Object.keys(q).length > 0 ? state.dataset.api.getDataPromise(q) : Promise.resolve({});
        return dataPromise.then(result => {
            dispatch(_setDotPlotData(updateDotPlotData(result.dotplot, dotPlotData, searchTokens)));
            const newSummary = result.summary || {};
            for (let key in newSummary) {
                globalFeatureSummary[key] = newSummary[key];
            }
            dispatch(setGlobalFeatureSummary(globalFeatureSummary));

            if (result.embedding) {
                for (let i = 0; i < result.embedding.length; i++) {
                    dispatch(handleEmbeddingResult(result.embedding[i], embeddingsUsed[i]));
                }
            }
            // embeddingData.sort((a, b) => {
            //     a = a.name.toLowerCase();
            //     b = b.name.toLowerCase();
            //     return a < b ? -1 : 1;
            // });

            if (state.chartOptions.activeEmbedding != null) { // when restoring view - put last so that it becomes active
                let index = -1;
                for (let i = 0; i < embeddingData.length; i++) {
                    if (state.chartOptions.activeEmbedding === getTraceKey(embeddingData[i])) {
                        index = i;
                        break;
                    }
                }
                if (index !== -1) {
                    let activeTrace = embeddingData[index];
                    embeddingData.splice(index, 1);
                    embeddingData.push(activeTrace);
                }
                state.chartOptions.activeEmbedding = null;
            }
            dispatch(setEmbeddingData(embeddingData.slice()));
            if (result.selection) {
                dispatch(handleSelectionResult(result.selection));
            } else {
                dispatch(_setSelectedDotPlotData(updateDotPlotData(null, selectedDotPlotData, searchTokens)));
            }
        }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve data. Please try again.');
            if (onError) {
                onError(err);
            }
        });

    };

}

// depends on global feature summary
function handleEmbeddingResult(embedding, embeddingDef) {
    return function (dispatch, getState) {
        const state = getState();

        // let embeddingResult = result.embeddingResult;

        let isSpatial = embeddingDef.spatial != null;
        let interpolator = state.interpolator;


        let embeddingData = state.embeddingData;
        const globalFeatureSummary = state.globalFeatureSummary;
        const obsCat = state.dataset.obsCat;
        const embeddingBins = embedding.bins;
        const embeddingValues = embedding.values;
        const coordinates = embedding.coordinates;


        let coordinateNames = [null, null, null]; // x, y, z
        for (let name in coordinates) {
            if (name.endsWith('_1')) {
                coordinateNames[0] = name;
            } else if (name.endsWith('_2')) {
                coordinateNames[1] = name;
            } else if (name.endsWith('_3')) {
                coordinateNames[2] = name;
            }
        }
        let x = coordinates[coordinateNames[0]];
        let y = coordinates[coordinateNames[1]];
        let z = coordinateNames[2] ? coordinates[coordinateNames[2]] : null;
        const is3d = z != null;

        // add new embedding values
        for (let name in embeddingValues) {
            let traceSummary = globalFeatureSummary[name];

            let values = embeddingValues[name];
            let purity = null;
            if (values.value !== undefined) {
                purity = values.purity;
                values = values.value;
            }

            let isCategorical = name !== '__count' && obsCat.indexOf(name) !== -1;
            let colorScale = null;

            if (!isCategorical) {
                // __count range is per embedding so we always recompute for now
                if (traceSummary == null || name === '__count') {

                    let min = Number.MAX_VALUE;
                    let max = -Number.MAX_VALUE;
                    for (let i = 0, n = values.length; i < n; i++) {
                        let value = values[i];
                        min = value < min ? value : min;
                        max = value > max ? value : max;
                    }
                    traceSummary = {min: min, max: max};
                    globalFeatureSummary[name] = traceSummary;
                }
                let domain = [traceSummary.min, traceSummary.max];
                if (traceSummary.customMin != null && !isNaN(traceSummary.customMin)) {
                    domain[0] = traceSummary.customMin;
                }
                if (traceSummary.customMax != null && !isNaN(traceSummary.customMax)) {
                    domain[1] = traceSummary.customMax;
                }
                colorScale = scaleSequential(interpolator.value).domain(domain);

            } else {
                let traceUniqueValues = traceSummary.categories;
                // if (traceUniqueValues.length === 1 && traceUniqueValues[0] === true) {
                //     traceUniqueValues = traceUniqueValues.concat([null]);
                //     traceSummary.categories = traceUniqueValues;
                //     traceSummary.counts = traceSummary.counts.concat([state.dataset.shape[0] - traceSummary.counts[0]]);
                // }

                let colors;
                if (traceUniqueValues.length <= 10) {
                    colors = schemeCategory10;
                } else if (traceUniqueValues.length <= 12) {
                    colors = schemePaired;
                } else if (traceUniqueValues.length <= 20) {
                    colors = CATEGORY_20B;
                } else {
                    colors = CATEGORY_20B.concat(CATEGORY_20C);
                }

                colorScale = scaleOrdinal(colors).domain(traceUniqueValues);
                colorScale.summary = traceSummary;
            }


            let chartData = {
                embedding: Object.assign({}, embeddingDef),
                name: name,
                x: x,
                y: y,
                z: is3d ? z : undefined,
                bins: embeddingBins,
                dimensions: is3d ? 3 : 2,
                npoints: x.length,
                date: new Date(),
                active: true,
                colorScale: colorScale,
                continuous: !isCategorical,
                isCategorical: isCategorical,
                values: values, // for color
                isImage: isSpatial
                // purity: purity,
                // text: values,
            };
            if (!isSpatial) {
                chartData.positions = getPositions(chartData);
            }
            updateTraceColors(chartData);

            if (isSpatial) {
                chartData.indices = !isCategorical ? indexSort(chartData.values, true) : randomSeq(chartData.values.length);
                chartData.isImage = true;
                const url = getState().dataset.api.getFileUrl(embeddingDef.spatial.image);
                let tileSource = new OpenSeadragon.ImageTileSource({
                    url: url,
                    buildPyramid: true,
                    crossOriginPolicy: "Anonymous"
                });
                chartData.tileSource = tileSource;
            }
            embeddingData.push(chartData);

        }

    };

}

function handleError(dispatch, err, message) {
    console.log(err);
    if (err.status === 401) {
        dispatch(setMessage('Your session has expired. Please login again.'));
        return dispatch(logout());
    }
    if (message == null) {
        message = err instanceof CustomError ? err.message : 'An unexpected error occurred. Please try again.';
    }
    dispatch(setMessage(new Error(message)));
}

export function listDatasets() {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        return getState().serverInfo.api.getDatasetsPromise()
            .then(choices => {
                dispatch(_setDatasetChoices(choices));
            })
            .finally(() => {
                dispatch(_setLoading(false));
            })
            .catch(err => {
                handleError(dispatch, err, 'Unable to retrieve datasets. Please try again.');
            });
    };
}



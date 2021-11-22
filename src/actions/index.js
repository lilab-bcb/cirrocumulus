import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10, schemePaired} from 'd3-scale-chromatic';
import {saveAs} from 'file-saver';
import {find, findIndex, groupBy, indexOf, isArray, isString} from 'lodash';
import OpenSeadragon from 'openseadragon';
import isPlainObject from 'react-redux/lib/utils/isPlainObject';
import CustomError from '../CustomError';
import {getPassingFilterIndices} from '../dataset_filter';
import {DirectAccessDataset} from '../DirectAccessDataset';
import {createCategoryToStats} from '../MetaEmbedding';

import {RestDataset} from '../RestDataset';
import {RestServerApi} from '../RestServerApi';
import {StaticServerApi} from '../StaticServerApi';

import {getPositions} from '../ThreeUtil';
import {
    addFeatureSetsToX,
    CATEGORY_20B,
    CATEGORY_20C,
    createColorScale,
    createEmbeddingDensity,
    FEATURE_TYPE,
    FEATURE_TYPE_MEASURES_EXCLUDE,
    getFeatureSets,
    getInterpolator,
    indexSort,
    randomSeq,
    SERVER_CAPABILITY_ADD_DATASET,
    SERVER_CAPABILITY_RENAME_CATEGORIES,
    summarizeDensity,
    TRACE_TYPE_IMAGE,
    TRACE_TYPE_META_IMAGE,
    TRACE_TYPE_SCATTER,
    updateTraceColors
} from '../util';
import {updateJob} from '../DotPlotJobResultsPanel';

export const API = process.env.REACT_APP_API_URL || 'api';
const authScopes = [
    'email'
    // 'profile',
    // 'https://www.googleapis.com/auth/userinfo.profile',
    // 'https://www.googleapis.com/auth/contacts.readonly',
    // 'https://www.googleapis.com/auth/devstorage.full_control',
];


export const DEFAULT_POINT_SIZE = 1;
export const DEFAULT_MARKER_OPACITY = 1;
export const DEFAULT_UNSELECTED_MARKER_OPACITY = 0.1;
export const DEFAULT_INTERPOLATORS = {};
DEFAULT_INTERPOLATORS[FEATURE_TYPE.X] = {name: 'Viridis', reversed: false, value: getInterpolator('Viridis')};
DEFAULT_INTERPOLATORS[FEATURE_TYPE.COUNT] = {name: 'Greys', reversed: false, value: getInterpolator('Greys')};
DEFAULT_INTERPOLATORS[FEATURE_TYPE.OBS] = {name: 'Inferno', reversed: false, value: getInterpolator('Inferno')};
DEFAULT_INTERPOLATORS[FEATURE_TYPE.MODULE] = {name: 'RdBu', reversed: true, value: getInterpolator('RdBu')};


export const DEFAULT_DISTRIBUTION_PLOT_INTERPOLATOR = 'Reds';
export const DEFAULT_DRAG_MODE = 'pan';

export const DEFAULT_SHOW_AXIS = true;
export const DEFAULT_SHOW_FOG = false;
export const DEFAULT_DARK_MODE = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)').matches : false;
export const DEFAULT_LABEL_FONT_SIZE = 14;
export const DEFAULT_LABEL_STROKE_WIDTH = 4;

export const SET_DRAG_DIVIDER = 'SET_DRAG_DIVIDER';
export const SET_WINDOW_SIZE = 'SET_WINDOW_SIZE';
export const SET_DRAWER_OPEN = 'SET_DRAWER_OPEN';
export const SET_EMBEDDING_LABELS = 'SET_EMBEDDING_LABELS';
export const SET_DISTRIBUTION_PLOT_INTERPOLATOR = 'SET_DISTRIBUTION_PLOT_INTERPOLATOR';
export const SET_CHART_OPTIONS = 'SET_CHART_OPTIONS';
export const SET_COMBINE_DATASET_FILTERS = 'SET_COMBINE_DATASET_FILTERS';
export const SET_DATASET_FILTERS = 'SET_DATASET_FILTERS'; // saved dataset filters
export const SET_DATASET_VIEWS = 'SET_DATASET_VIEWS'; // saved dataset views


export const SET_LEGEND_SCROLL_POSITION = 'SET_LEGEND_SCROLL_POSITION';
export const SET_ACTIVE_FEATURE = 'SET_ACTIVE_FEATURE';
export const SET_CHART_SIZE = 'SET_CHART_SIZE';
export const SET_SERVER_INFO = "SET_SERVER_INFO";
export const SET_DATASET_FILTER = 'SET_DATASET_FILTER';
export const ADD_DATASET = 'ADD_DATASET';
export const DELETE_DATASET = 'DELETE_DATASET';
export const UPDATE_DATASET = 'UPDATE_DATASET';
export const SET_GLOBAL_FEATURE_SUMMARY = 'SET_GLOBAL_FEATURE_SUMMARY';
export const SET_SAVED_DATASET_STATE = 'SET_SAVED_DATASET_STATE';

export const SET_DOMAIN = 'SET_DOMAIN';
export const UPDATE_CATEGORICAL_COLOR = 'UPDATE_CATEGORICAL_COLOR';
export const SET_CATEGORICAL_NAME = 'SET_CATEGORICAL_NAME';
export const UPDATE_CATEGORICAL_NAME = 'UPDATE_CATEGORICAL_NAME';
export const SET_MARKER_SIZE = 'SET_MARKER_SIZE';
export const SET_MARKER_OPACITY = 'SET_MARKER_OPACITY';


export const SET_UNSELECTED_MARKER_OPACITY = 'SET_UNSELECTED_MARKER_OPACITY';

export const SET_SELECTION = 'SET_SELECTION';
export const SET_FEATURE_SUMMARY = 'SET_FEATURE_SUMMARY';

export const SET_SEARCH_TOKENS = 'SET_SEARCH_TOKENS';

export const SET_SELECTED_EMBEDDING = 'SET_SELECTED_EMBEDDING';
export const SET_MESSAGE = 'SET_MESSAGE';
export const SET_INTERPOLATOR = 'SET_INTERPOLATOR';
export const SET_POINT_SIZE = 'SET_POINT_SIZE';
export const SET_UNSELECTED_POINT_SIZE = 'SET_UNSELECTED_POINT_SIZE';

export const SET_EMAIL = 'SET_EMAIL';
export const SET_USER = 'SET_USER';
export const SET_DATASET = 'SET_DATASET';
export const SET_MARKERS = 'SET_MARKERS';
export const SET_DIALOG = 'SET_DIALOG';

export const OPEN_DATASET_DIALOG = 'OPEN_DATASET_DIALOG';
export const EDIT_DATASET_DIALOG = 'EDIT_DATASET_DIALOG';
export const IMPORT_DATASET_DIALOG = 'IMPORT_DATASET_DIALOG';
export const SAVE_DATASET_FILTER_DIALOG = 'SAVE_DATASET_FILTER_DIALOG';
export const SAVE_FEATURE_SET_DIALOG = 'SAVE_FEATURE_SET_DIALOG';
export const HELP_DIALOG = 'HELP_DIALOG';
export const DELETE_DATASET_DIALOG = 'DELETE_DATASET_DIALOG';

export const SET_DATASET_CHOICES = 'SET_DATASET_CHOICES';
export const RESTORE_VIEW = 'RESTORE_VIEW';

export const SET_DISTRIBUTION_DATA = 'SET_DISTRIBUTION_DATA';
export const SET_SELECTED_DISTRIBUTION_DATA = 'SET_SELECTED_DISTRIBUTION_DATA';

export const SET_DISTRIBUTION_PLOT_OPTIONS = 'SET_DISTRIBUTION_PLOT_OPTIONS';
export const SET_EMBEDDING_DATA = 'SET_EMBEDDING_DATA';

export const ADD_TASK = 'ADD_TASK';
export const REMOVE_TASK = 'REMOVE_TASK';

export const SET_TAB = 'SET_TAB';

export const SET_LOADING_APP = 'LOADING_APP';

export const SET_JOB_RESULTS = 'SET_JOB_RESULTS';
export const SET_JOB_RESULT = 'SET_JOB_RESULT';


export function getEmbeddingKey(embedding, includeDensity = true) {
    let key = embedding.name;
    if (embedding.mode != null && includeDensity) {
        key += '_' + embedding.mode;
    }
    return key;
}


export function getTraceKey(traceInfo) {
    return traceInfo.name + '_' + getEmbeddingKey(traceInfo.embedding);
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
        const approximateColdBootTime = 10;

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
        if (process.env.REACT_APP_STATIC === 'true') {
            const serverInfo = {clientId: '', capabilities: new Set()};
            serverInfo.api = new StaticServerApi();
            dispatch(setServerInfo(serverInfo));
            dispatch(_setLoadingApp({loading: false}));
            dispatch(listDatasets()).then(() => {
                dispatch(_loadSavedView());
            });
            return Promise.resolve();
        }
        let urlOptions = window.location.search.substring(1);
        // can be run in server mode by passing url not stored in database
        let externalDatasetId = null;
        if (urlOptions.length > 0) {
            let pairs = urlOptions.split('&');
            pairs.forEach(pair => {
                let keyVal = pair.split('=');
                if (keyVal.length === 2 && (keyVal[0] === 'url' || keyVal[0] === 'id')) { // load a dataset by URL or id
                    externalDatasetId = keyVal[1];
                }
            });
        }
        const getServerInfo = fetch(API + '/server').then(result => result.json());
        return getServerInfo.then(serverInfo => {
            serverInfo.api = new RestServerApi();
            if (serverInfo.clientId == null) {
                serverInfo.clientId = '';
            }
            const capabilities = new Set();
            for (let key in serverInfo.capabilities) {
                if (serverInfo.capabilities[key]) {
                    capabilities.add(key);
                }
            }
            serverInfo.capabilities = capabilities;
            dispatch(setServerInfo(serverInfo));
            if (serverInfo.clientId === '' || serverInfo.clientId == null) { // no login required
                dispatch(_setLoadingApp({loading: false}));
                dispatch(_setEmail(capabilities.has(SERVER_CAPABILITY_ADD_DATASET) ? '' : null));
                if (capabilities.has(SERVER_CAPABILITY_ADD_DATASET)) {
                    dispatch(setUser({importer: true}));
                }
                if (externalDatasetId == null) {
                    dispatch(listDatasets()).then(() => {
                        dispatch(_loadSavedView());
                    });
                } else {
                    dispatch(_setDatasetChoices([{
                        id: externalDatasetId,
                        url: externalDatasetId,
                        name: externalDatasetId,
                        description: ''
                    }]));
                    dispatch(_loadSavedView());
                }

            } else {
                console.log((new Date().getTime() - startTime) / 1000 + " startup time");
                let script = document.createElement('script');
                script.type = 'text/javascript';
                script.src = 'https://apis.google.com/js/api.js';
                script.onload = (e) => {
                    window.gapi.load('client:auth2', () => {
                        window.gapi.client.init({
                            clientId: serverInfo.clientId,
                            scope: authScopes.join(' ')
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


export function openView(id, loadDataset = false) {
    return function (dispatch, getState) {
        const task = {name: 'Open view'};
        dispatch(addTask(task));
        getState().serverInfo.api.getViewPromise(id)
            .then(result => {
                const savedView = isString(result.value) ? JSON.parse(result.value) : result.value;
                savedView.dataset = result.dataset_id;
                dispatch(restoreSavedView(savedView));
            }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve link. Please try again.');
        });
    };
}


/**
 *
 * @param payload Object with name, notes
 * @returns {(function(*=, *): void)|*}
 */
export function saveView(payload) {
    return function (dispatch, getState) {
        const state = getState();
        const value = getDatasetStateJson(state);
        payload.value = value;
        delete value['dataset'];
        payload = Object.assign({ds_id: state.dataset.id}, payload);
        const task = {name: 'Save view'};
        dispatch(addTask(task));
        getState().serverInfo.api.upsertViewPromise(payload, false)
            .then(result => {
                payload = Object.assign(result, payload);
                const array = getState().datasetViews;
                array.push(payload);
                dispatch(setDatasetViews(array.slice()));
                dispatch(setMessage('Link saved'));
            }).finally(() => {
            dispatch(removeTask(task));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to save link. Please try again.');
        });
    };
}

export function deleteView(id) {
    return function (dispatch, getState) {
        const task = {name: 'Delete view'};
        dispatch(addTask(task));
        getState().serverInfo.api.deleteViewPromise(id, getState().dataset.id)
            .then(() => {
                let array = getState().datasetViews;
                for (let i = 0; i < array.length; i++) {
                    if (array[i].id === id) {
                        array.splice(i, 1);
                        break;
                    }
                }
                dispatch(setDatasetFilters(array.slice()));
                dispatch(setMessage('Link deleted'));
            }).finally(() => {
            dispatch(removeTask(task));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete link. Please try again.');
        });
    };
}


export function submitJob(jobData) {
    return function (dispatch, getState) {
        let jobId;
        let timeout = 5 * 1000; // TODO
        function getJobStatus() {
            getState().dataset.api.getJobStatus(jobId)
                .then(result => {
                    const jobResult = find(getState().jobResults, item => item.id === jobId);
                    if (jobResult == null) { // job was deleted
                        return;
                    }
                    const statusUpdated = jobResult.status !== result.status;
                    jobResult.status = result.status;
                    let fetchJobStatus = true;
                    if (result.status === 'complete') {
                        dispatch(setMessage(jobData.name + ': complete'));
                        dispatch(setJobResults(getState().jobResults.slice()));
                        fetchJobStatus = false;
                    } else if (result.status === 'error') {
                        handleError(dispatch, new CustomError('Unable to complete job. Please try again.'));
                        fetchJobStatus = false;
                    } else if (statusUpdated) { // force repaint
                        dispatch(setMessage(jobData.name + ': ' + result.status));
                        dispatch(setJobResults(getState().jobResults.slice()));
                    }
                    if (fetchJobStatus) {
                        window.setTimeout(getJobStatus, timeout);
                    }
                }).catch(err => {
                handleError(dispatch, err, 'Unable to get job. Please try again.');
            });
        }

        jobData.id = getState().dataset.id;
        if (jobData.params.filter == null) {
            jobData.params.filter = getFilterJson(getState());
        }
        getState().serverInfo.api.submitJob(jobData)
            .then(result => {
                dispatch(setMessage('Job submitted'));
                jobId = result.id;
                jobData.email = getState().email;
                jobData.id = jobId;
                getState().jobResults.push(jobData);
                dispatch(setJobResults(getState().jobResults.slice()));
                window.setTimeout(getJobStatus, timeout);
            }).finally(() => {

        }).catch(err => {
            handleError(dispatch, err, 'Unable to submit job. Please try again.');
        });
    };
}

export function deleteFeatureSet(id) {
    return function (dispatch, getState) {
        const task = {name: 'Delete set'};
        dispatch(addTask(task));
        getState().serverInfo.api.deleteFeatureSet(id, getState().dataset.id)
            .then(result => {
                let markers = getState().markers;
                let found = false;
                for (let i = 0; i < markers.length; i++) {
                    if (markers[i].id === id) {
                        markers.splice(i, 1);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    console.log('Unable to find feature set id ' + id);
                }
                dispatch(setMarkers(markers.slice()));
                dispatch(setMessage('Set deleted'));
            }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete set. Please try again.');
        });
    };
}

export function saveFeatureSet(payload) {
    return function (dispatch, getState) {
        const state = getState();
        const features = state.searchTokens.filter(item => item.type === FEATURE_TYPE.X).map(item => item.id);
        const requestBody = {
            ds_id: state.dataset.id,
            name: payload.name,
            features: features,
            category: payload.category
        };
        const task = {name: 'Save feature set'};
        dispatch(addTask(task));
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
            dispatch(removeTask(task));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to save set. Please try again.');
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
                    if (Array.isArray(datasetFilter[key])) {
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

function getDatasetFilterDependencies(datasetFilter) {
    let features = new Set();
    let basis = new Set();
    for (let key in datasetFilter) {
        // basis, path for brush filter
        const filterObject = datasetFilter[key];
        if (Array.isArray(filterObject)) { // brush filter

        } else if (filterObject.operation === 'in') {
            features.add(key);
        } else {
            if (filterObject.operation !== '' && !isNaN(filterObject.value) && filterObject.value != null) {
                features.add(key);
            }
        }
    }
    return {features: features, basis: basis};
}

export function getDatasetFilterNames(datasetFilter) {
    const names = [];
    let isBrushing = false;
    for (let key in datasetFilter) {
        const value = datasetFilter[key];
        if (Array.isArray(value)) {
            isBrushing = true;
        } else if (value.operation === 'in') {
            names.push(key);
        } else if (value.operation != null && value.operation !== '' && !isNaN(value.value) && value.value != null) {
            names.push(key);
        }
    }
    if (isBrushing) {
        names.push('selection');
    }
    return names;
}

export function getDatasetFilterArray(datasetFilter) {
    let filters = [];
    const brushIndices = new Set();
    // brush filters are combined with OR
    for (let key in datasetFilter) {
        // basis, path for brush filter
        const value = datasetFilter[key];
        if (Array.isArray(value)) {
            value.forEach(brushFilter => {
                brushFilter.indices.forEach(index => brushIndices.add(index));
            });
        }
    }
    if (brushIndices.size > 0) {
        filters.push(['__index', 'in', brushIndices]);
    }
    for (let key in datasetFilter) {
        const value = datasetFilter[key];
        let f = null;
        if (Array.isArray(value)) {
            continue;
        }
        if (value.operation === 'in') {
            f = [key, value.operation, value.value];
        } else if (value.operation != null && value.operation !== '' && !isNaN(value.value) && value.value != null) {
            f = [key, value.operation, value.value];
        }
        if (f != null) {
            filters.push(f);
        }
    }
    return filters;
}


function getFilterJson(state) {
    return datasetFilterToJson(state.dataset, state.datasetFilter, state.combineDatasetFilters);
}

export function datasetFilterToJson(dataset, datasetFilter, combineDatasetFilters) {
    let filters = getDatasetFilterArray(datasetFilter);
    if (filters.length > 0) {
        const obs = dataset.obs;
        const obsCat = dataset.obsCat;
        for (let i = 0; i < filters.length; i++) {
            // add obs/ prefix
            const filter = filters[i];
            if (filter[0] === '__index') {
                filter[2] = Array.from(filter[2]); // convert Set to array
            } else if (obsCat.indexOf(filter[0]) !== -1 || obs.indexOf(filter[0]) !== -1) {
                filter[0] = 'obs/' + filter[0];
            }
        }
        return {filters: filters, combine: combineDatasetFilters};
    }
}


export function downloadSelectedIds() {
    return function (dispatch, getState) {
        const task = {name: 'Download ids'};
        dispatch(addTask(task));
        const state = getState();
        let filter = getFilterJson(state, true);
        state.dataset.api.getSelectedIdsPromise({
            filter: filter
        }, state.cachedData).then(result => {
            const blob = new Blob([result.ids.join('\n')], {type: "text/plain;charset=utf-8"});
            saveAs(blob, "selection.txt");
        }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err);
        });
    };
}

export function exportDatasetFilters() {
    return function (dispatch, getState) {
        const task = {name: 'Export Filters'};
        dispatch(addTask(task));
        getState().serverInfo.api.exportDatasetFiltersPromise(getState().dataset.id).then(result => {
            if (result == null) {
                handleError(dispatch, 'Unable to export filters');
                return;
            }
            const blob = new Blob([result], {type: "text/plain;charset=utf-8"});
            saveAs(blob, "filters.csv");
        }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err);
        });
    };
}


export function setDistributionPlotOptions(payload) {
    return {type: SET_DISTRIBUTION_PLOT_OPTIONS, payload: payload};
}

export function setChartOptions(payload) {
    return {type: SET_CHART_OPTIONS, payload: payload};
}

export function setDrawerOpen(payload) {
    return {type: SET_DRAWER_OPEN, payload: payload};
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


export function setWindowSize(payload) {
    return {type: SET_WINDOW_SIZE, payload: payload};
}

export function setDragDivider(payload) {
    return {type: SET_DRAG_DIVIDER, payload: payload};
}

export function _setActiveFeature(payload) {
    return {type: SET_ACTIVE_FEATURE, payload: payload};
}


export function setActiveFeature(payload) {
    return {type: SET_ACTIVE_FEATURE, payload: payload};
}

export function setLegendScrollPosition(payload) {
    return {type: SET_LEGEND_SCROLL_POSITION, payload: payload};
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
        const task = {name: 'Update Filter'};
        dispatch(addTask(task));
        // whenever filter is updated, we need to get selection statistics
        const state = getState();
        const searchTokens = state.searchTokens;
        let filter = getFilterJson(state);
        const groupedSearchTokens = groupBy(searchTokens, 'type');
        addFeatureSetsToX(getFeatureSets(state.markers, groupedSearchTokens[FEATURE_TYPE.FEATURE_SET] || []), (searchTokens[FEATURE_TYPE.X] || []).map(item => item.id));

        const measures = [];
        for (const key in groupedSearchTokens) {
            if (FEATURE_TYPE_MEASURES_EXCLUDE.indexOf(key) === -1) {
                const prefix = key === FEATURE_TYPE.X ? '' : key + '/';
                groupedSearchTokens[key].forEach(item => measures.push(prefix + item.id));
            }
        }

        let q = {
            selection: {
                measures: measures,
                dimensions: (groupedSearchTokens[FEATURE_TYPE.OBS_CAT] || []).map(item => item.id)
            }
        };

        if (filter) {
            q.selection.filter = filter;
        }

        if (filter == null) {
            if (state.selection != null) { // reset
                dispatch(setSelection(null));
            }
            dispatch(setSelectedDistributionData([]));
            dispatch(setFeatureSummary({}));
            dispatch(removeTask(task));
            return;
        }

        const cachedData = state.cachedData;
        getState().dataset.api.getDataPromise(q, cachedData).then(result => {
            dispatch(handleSelectionResult(result.selection, true));
        }).catch(err => {
            handleError(dispatch, err);
        }).finally(() => dispatch(removeTask(task)));
    };
}


export function handleBrushFilterUpdated(payload) {
    return function (dispatch, getState) {
        const name = payload.name; // full basis name
        const value = payload.value;  // value has basis and indices
        const clear = payload.clear;
        const indices = payload.value != null ? payload.value.indices : null;
        let datasetFilter = getState().datasetFilter;

        function clearFilters() {
            if (clear) {
                for (let key in datasetFilter) {
                    const value = datasetFilter[key];
                    if (Array.isArray(value)) {
                        delete datasetFilter[key];
                    }
                }
            }
        }

        let update = true;
        if (value == null || indices.size === 0) { // remove filter
            update = datasetFilter[name] != null;
            delete datasetFilter[name];
            clearFilters();
        } else {
            let priorFilters = datasetFilter[name];
            let isToggleRegion = false;
            if (priorFilters != null) {
                const priorIndex = value.id != null ? findIndex(priorFilters, f => f.id === value.id) : -1;
                if (priorIndex !== -1) { // toggle region
                    isToggleRegion = true;
                    priorFilters.splice(priorIndex, 1);
                }
            }
            if (!isToggleRegion) {
                clearFilters();
                priorFilters = datasetFilter[name];
                if (priorFilters == null) {
                    priorFilters = [];
                    datasetFilter[name] = priorFilters;
                }
                priorFilters.push({basis: value.basis, indices: indices, id: value.id});
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

export function handleDimensionFilterUpdated(payload) {
    return function (dispatch, getState) {
        let name = payload.name;
        let newValue = payload.value;
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
        let categoricalFilter = datasetFilter[name];
        if (categoricalFilter == null) {
            categoricalFilter = {operation: 'in', value: []};
            datasetFilter[name] = categoricalFilter;
        }

        if (shiftKey && categoricalFilter.value.length > 0) {
            // add last click to current
            let lastIndex = categories.indexOf(categoricalFilter.value[categoricalFilter.value.length - 1]);
            let currentIndex = categories.indexOf(newValue);
            // put clicked category at end of array
            if (currentIndex > lastIndex) {
                for (let i = lastIndex; i <= currentIndex; i++) {
                    let index = categoricalFilter.value.indexOf(newValue);
                    if (index !== -1) {
                        categoricalFilter.value.splice(index, 1);
                    }
                    categoricalFilter.value.push(categories[i]);
                }
            } else {
                for (let i = lastIndex; i >= currentIndex; i--) {
                    let index = categoricalFilter.value.indexOf(newValue);
                    if (index !== -1) {
                        categoricalFilter.value.splice(index, 1);
                    }
                    categoricalFilter.value.push(categories[i]);
                }
            }
        } else {
            let selectedIndex = categoricalFilter.value.indexOf(newValue);
            if (!metaKey) { // clear and toggle current
                categoricalFilter.value = [];
            }
            if (selectedIndex !== -1) { // exists, remove
                categoricalFilter.value.splice(selectedIndex, 1);
                if (categoricalFilter.value.length === 0) {
                    delete datasetFilter[name];
                }
            } else {
                categoricalFilter.value.push(newValue);
            }
        }
        dispatch(setDatasetFilter(datasetFilter));
        dispatch(handleFilterUpdated());
    };
}


export function handleDomainChange(payload) {
    return {type: SET_DOMAIN, payload: payload};
}

function _handleCategoricalNameChange(payload) {
    return {type: SET_CATEGORICAL_NAME, payload: payload};
}

function handleUpdateCategoricalName(payload) {
    return {type: UPDATE_CATEGORICAL_NAME, payload: payload};
}

export function _handleColorChange(payload) {
    return {type: UPDATE_CATEGORICAL_COLOR, payload: payload};
}

export function handleColorChange(payload) {
    return function (dispatch, getState) {
        const value = Object.assign({id: getState().dataset.id}, payload);
        if (getState().serverInfo.capabilities.has(SERVER_CAPABILITY_RENAME_CATEGORIES)) {
            getState().serverInfo.api.setCategoryNamePromise(value).then(() => {
                dispatch(_handleColorChange(payload));
            });
        } else { // save, but do not persist
            dispatch(_handleColorChange(payload));
        }
    };
}

export function handleCategoricalNameChange(payload) {
    return function (dispatch, getState) {
        const value = Object.assign({id: getState().dataset.id}, payload);
        if (getState().serverInfo.capabilities.has(SERVER_CAPABILITY_RENAME_CATEGORIES)) {
            getState().serverInfo.api.setCategoryNamePromise(value).then(() => {
                dispatch(handleUpdateCategoricalName(payload));
            });
        } else { // save, but do not persist
            dispatch(handleUpdateCategoricalName(payload));
        }
    };
}


export function restoreView(payload) {
    return {type: RESTORE_VIEW, payload: payload};
}


export function setPointSize(payload) {
    return {type: SET_POINT_SIZE, payload: payload};
}

export function setUnselectedPointSize(payload) {
    return {type: SET_UNSELECTED_POINT_SIZE, payload: payload};
}

export function setMarkerSize(payload) {
    return {type: SET_MARKER_SIZE, payload: payload};
}


export function setUser(payload) {
    return {type: SET_USER, payload: payload};
}


export function setServerInfo(payload) {
    return {type: SET_SERVER_INFO, payload: payload};
}

function getDefaultDatasetView(dataset) {
    const embeddingNames = dataset.embeddings.map(e => e.name);
    let selectedEmbedding = null;
    let obsCat = null;
    if (embeddingNames.length > 0) {
        let embeddingPriorities = ['tissue_hires', 'fle', 'umap', 'tsne'];
        let embeddingName = null;
        for (let priorityIndex = 0; priorityIndex < embeddingPriorities.length && embeddingName == null; priorityIndex++) {
            for (let i = 0; i < embeddingNames.length; i++) {
                if (embeddingNames[i].toLowerCase().indexOf(embeddingPriorities[priorityIndex]) !== -1) {
                    embeddingName = embeddingNames[i];
                    break;
                }
            }
        }

        if (embeddingName == null) {
            embeddingName = embeddingNames[0];
        }
        selectedEmbedding = dataset.embeddings[dataset.embeddings.map(e => e.name).indexOf(embeddingName)];
    }
    if (dataset.markers != null && dataset.markers.length > 0) {
        let category = dataset.markers[0].category;
        const suffix = ' (rank_genes_';
        let index = category.indexOf(suffix);
        if (index !== -1) {
            category = category.substring(0, index);
        }
        if (dataset.obsCat.indexOf(category) !== -1) {
            obsCat = category;
        }
    }
    if (obsCat == null) {
        let catPriorities = ['anno', 'cell_type', 'celltype', 'leiden', 'louvain', 'seurat_cluster', 'cluster'];
        for (let priorityIndex = 0; priorityIndex < catPriorities.length && obsCat == null; priorityIndex++) {
            for (let i = 0; i < dataset.obsCat.length; i++) {
                if (dataset.obsCat[i].toLowerCase().indexOf(catPriorities[priorityIndex]) !== -1) {
                    obsCat = dataset.obsCat[i];
                    break;
                }
            }
        }
    }

    return {selectedEmbedding, obsCat};
}

function loadDefaultDatasetView() {
    return function (dispatch, getState) {
        const dataset = getState().dataset;
        const {selectedEmbedding, obsCat} = getDefaultDatasetView(dataset);
        if (obsCat != null) {
            dispatch(setSearchTokens([{id: obsCat, type: FEATURE_TYPE.OBS_CAT}]));
        }
        if (selectedEmbedding != null) {
            dispatch(setSelectedEmbedding([selectedEmbedding]));
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
        if (savedView.interpolator != null) {
            for (const key in savedView.interpolator) {
                savedView.interpolator[key].value = getInterpolator(savedView.interpolator[key].name);
            }
        }

        if (savedView.datasetFilter == null) {
            savedView.datasetFilter = {};

        } else {
            for (let key in savedView.datasetFilter) {
                let value = savedView.datasetFilter[key];
                if (value.operation) {
                    value.uiValue = value.value;
                }
            }
        }
        const task = {name: 'Restore view'};
        dispatch(addTask(task));
        let datasetPromise;
        if (savedView.dataset != null) {
            datasetPromise = dispatch(setDataset(savedView.dataset, false, false));
        } else {
            datasetPromise = Promise.resolve();
        }
        datasetPromise
            .then(() => {
                let dataset = getState().dataset;
                if (savedView.embeddings && savedView.embeddings.length > 0) {
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
                        if (savedView.chartOptions == null) {
                            savedView.chartOptions = {};
                        }
                        savedView.chartOptions.camera = savedView.camera;
                    }
                } else {
                    const {selectedEmbedding, obsCat} = getDefaultDatasetView(dataset);
                    if (selectedEmbedding) {
                        savedView.embeddings = [selectedEmbedding];
                        if ((savedView.q == null || savedView.q.length === 0) && obsCat != null) {
                            savedView.q = [{id: obsCat, type: FEATURE_TYPE.OBS_CAT}];
                            savedView.activeFeature = {
                                name: obsCat,
                                type: FEATURE_TYPE.OBS_CAT,
                                embeddingKey: obsCat + '_' + getEmbeddingKey(selectedEmbedding)
                            };
                        }
                    }

                }
            })
            .then(() => dispatch(setDatasetFilter(savedView.datasetFilter)))
            .then(() => dispatch(restoreView(savedView)))
            .then(() => dispatch(_updateCharts()))
            .then(() => dispatch(handleFilterUpdated()))
            .then(() => {
                if (savedView.distributionPlotOptions != null) {
                    dispatch(setDistributionPlotOptions(savedView.distributionPlotOptions));
                }
                let activeFeature = savedView.activeFeature;

                if (activeFeature == null && savedView.embeddings && savedView.embeddings.length > 0 && savedView.q && savedView.q.length > 0) { // pick the 1st search token and 1st embedding
                    activeFeature = {
                        name: savedView.q[0].id,
                        type: savedView.q[0].type,
                        embeddingKey: savedView.q[0].id + '_' + getEmbeddingKey(savedView.embeddings[0])
                    };
                }
                if (activeFeature != null) {
                    dispatch(setActiveFeature(activeFeature));
                }
            })
            .then(() => {
                if (savedView.jobId != null) {
                    dispatch(setJobResultId(savedView.jobId));
                }
            })
            .finally(() => dispatch(removeTask(task)))
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
        // #q=
        let q = window.location.hash.substring(3);
        if (q.length > 0) {
            try {
                savedView = JSON.parse(window.decodeURIComponent(q));
            } catch (err) {
                return dispatch(setMessage('Unable to restore view.'));
            }
        }
        if (savedView.link != null) {
            dispatch(openView(savedView.link, true));
        } else if (savedView.dataset != null) {
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
        let existingDataset = payload.dataset;
        const isEdit = existingDataset != null;
        const formData = {};
        if (existingDataset != null) {
            const existingReaders = new Set(existingDataset.readers);
            if (payload.readers != null) {
                let setsEqual = existingReaders.size === payload.readers.length;
                if (setsEqual) {
                    for (let i = 0; i < payload.readers.length; i++) {
                        if (!existingReaders.has(payload.readers[i])) {
                            setsEqual = false;
                            break;
                        }
                    }
                }
                if (!setsEqual) {
                    formData.readers = payload.readers;
                }
            }
        }
        if (existingDataset == null) {
            existingDataset = {};
        }
        const blacklist = new Set(['readers', 'dataset']);
        for (let key in payload) {
            if (!blacklist.has(key)) {
                const value = payload[key];
                if (value != null && value !== existingDataset[key]) {
                    formData[key] = value;
                }
            }
        }


        if (Object.keys(formData).length === 0) {
            dispatch(setDialog(null));
            return;
        }
        if (isEdit) {
            formData.id = payload.dataset.id;
        }
        const task = {name: 'Save Dataset'};
        dispatch(addTask(task));
        const request = getState().serverInfo.api.upsertDatasetPromise(formData);
        request.upload.addEventListener('progress', function (e) {
            if (formData['file'] != null) {
                let percent = (e.loaded / e.total) * 100;
                dispatch(setMessage('Percent ' + percent));
            }
        });
        request.addEventListener('load', function (e) {
            dispatch(removeTask(task));
            dispatch(setDialog(null));
            const status = request.status;
            if (status != 200) {
                return dispatch(setMessage('Unable to save dataset. Please try again.'));
            }

            const resp = JSON.parse(request.response);
            delete formData['file'];
            existingDataset = Object.assign({}, existingDataset, formData, resp);
            existingDataset.owner = true;

            if (isEdit) {
                dispatch(updateDataset(existingDataset));
                dispatch(setMessage('Dataset updated'));
            } else {
                dispatch(_addDataset(existingDataset));
                dispatch(setMessage('Dataset added'));
            }
        });
    };
}

export function deleteDataset(payload) {
    return function (dispatch, getState) {
        const task = {name: 'Delete Dataset'};
        dispatch(addTask(task));
        getState().serverInfo.api.deleteDatasetPromise(payload.dataset.id).then(() => {
            dispatch(_setDataset(null));
            dispatch(_deleteDataset({id: payload.dataset.id}));
            dispatch(setDialog(null));
            dispatch(setMessage('Dataset deleted'));
        }).finally(() => {
            dispatch(removeTask(task));
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

export function setDistributionPlotInterpolator(payload) {
    return {type: SET_DISTRIBUTION_PLOT_INTERPOLATOR, payload: payload};
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

function addTask(payload) {
    return {type: ADD_TASK, payload: payload};
}

function removeTask(payload) {
    return {type: REMOVE_TASK, payload: payload};
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

export function setJobResults(payload) {
    return {type: SET_JOB_RESULTS, payload: payload};
}

function _setJobResultId(payload) {
    return {type: SET_JOB_RESULT, payload: payload};
}

export function deleteJobResult(payload) {
    return function (dispatch, getState) {
        const task = {name: 'Delete Job'};
        dispatch(addTask(task));
        getState().dataset.api.deleteJob(payload).then(() => {
            let jobResults = getState().jobResults;
            const index = findIndex(jobResults, item => item.id === payload);
            jobResults.splice(index, 1);
            if (getState().jobResult === payload) {
                dispatch(_setJobResultId(null));
            }
            dispatch(setJobResults(jobResults.slice()));
        }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete result. Please try again.');
        });
    };
}

export function setJobResultId(jobId) {
    return function (dispatch, getState) {
        const existingJobResult = find(getState().jobResults, item => item.id === jobId);
        if (existingJobResult.data != null) { // data already loaded
            if (existingJobResult.type === 'de') {
                updateJob(existingJobResult);
            }
            return dispatch(_setJobResultId(existingJobResult.id));
        }
        const task = {name: 'Open Job'};
        dispatch(addTask(task));

        let promises = [];
        promises.push(getState().dataset.api.getJobParams(jobId));
        promises.push(getState().dataset.api.getJob(jobId));
        Promise.all(promises).then((values) => {
            const params = values[0];
            let result = values[1];
            const jobResults = getState().jobResults;
            const index = findIndex(jobResults, item => item.id === jobId);
            if (index === -1) { // job was deleted while fetching result?
                console.log('Unable to find job');
            } else {
                if (result.data == null) {
                    result = {data: result};
                }
                const jobResult = Object.assign({}, params, jobResults[index], result);
                if (jobResult.type === 'de') {
                    updateJob(jobResult);
                }
                let promise = Promise.resolve({});
                if (jobResult.params && jobResult.params.obs) {
                    const cachedData = getState().cachedData;
                    const q = {values: {dimensions: []}};
                    jobResult.params.obs.forEach(field => {
                        if (cachedData[field] == null) {
                            q.values.dimensions.push(field);
                        }
                    });
                    if (q.values.dimensions.length > 0) {
                        promise = getState().dataset.api.getDataPromise(q, cachedData);
                    }
                }
                jobResults[index] = jobResult;
                promise.then(() => {
                    dispatch(_setJobResultId(jobResult.id));
                });
            }
        }).finally(() => {
            dispatch(removeTask(task));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve result. Please try again.');
        });

    };
}

export function setMarkerOpacity(payload) {
    return {type: SET_MARKER_OPACITY, payload: payload};
}

export function setUnselectedMarkerOpacity(payload) {
    return {type: SET_UNSELECTED_MARKER_OPACITY, payload: payload};
}

export function toggleEmbeddingLabel(payload) {
    return function (dispatch, getState) {
        const embeddingLabels = getState().embeddingLabels;
        const index = embeddingLabels.indexOf(payload);
        if (index === -1) {
            embeddingLabels.push(payload);
        } else {
            embeddingLabels.splice(index, 1);
        }
        return dispatch({type: SET_EMBEDDING_LABELS, payload: embeddingLabels.slice()});
    };

}

function setDistributionData(payload) {
    return {type: SET_DISTRIBUTION_DATA, payload: payload};
}

function setSelectedDistributionData(payload) {
    return {type: SET_SELECTED_DISTRIBUTION_DATA, payload: payload};
}


function _setEmail(payload) {
    return {type: SET_EMAIL, payload: payload};
}


export function setSearchTokens(tokens, updateActiveFeature = true) {
    return function (dispatch, getState) {
        dispatch({type: SET_SEARCH_TOKENS, payload: tokens});
        dispatch(_updateCharts(null, updateActiveFeature));
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


export function setSavedDatasetState(payload) {
    return {type: SET_SAVED_DATASET_STATE, payload: payload};
}

function setDatasetFilters(payload) {
    return {type: SET_DATASET_FILTERS, payload: payload};
}

function setDatasetViews(payload) {
    return {type: SET_DATASET_VIEWS, payload: payload};
}


export function setDataset(id, loadDefaultView = true, setLoading = true) {
    return function (dispatch, getState) {
        let savedDatasetState = getState().savedDatasetState[id];
        const datasetChoices = getState().datasetChoices;
        let selectedChoice = null; // has id, owner, name
        for (let i = 0; i < datasetChoices.length; i++) {
            if (datasetChoices[i].id === id) {
                selectedChoice = datasetChoices[i];
                break;
            }
        }
        if (selectedChoice == null) { // search by name
            for (let i = 0; i < datasetChoices.length; i++) {
                if (datasetChoices[i].name === id) {
                    selectedChoice = datasetChoices[i];
                    break;
                }
            }
        }
        if (selectedChoice == null) { //
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

        let categoryNameResults;
        let datasetViews = [];
        let newDataset;
        let jobResults = [];

        function onPromisesComplete() {
            newDataset = Object.assign({}, dataset, newDataset);
            newDataset.api = dataset.api;
            newDataset.id = id;
            dispatch(_setDataset(newDataset));

            if (newDataset.results && newDataset.results.length > 0) {
                jobResults = jobResults.concat(newDataset.results);
            }
            dispatch(setJobResults(jobResults));
            // if (jobResults.length === 1) {
            //     dispatch(setJobResult(jobResults[0].id));
            // }

            if (categoryNameResults != null) {
                dispatch(_handleCategoricalNameChange(categoryNameResults));
            }
            dispatch(setDatasetViews(datasetViews));
            // dispatch(setDatasetFilters(datasetFilters));
            if (savedDatasetState) {
                dispatch(restoreSavedView(savedDatasetState));
            } else if (loadDefaultView) {
                dispatch(loadDefaultDatasetView());
            }
        }

        const task = setLoading ? {name: 'Set dataset'} : null;
        if (task) {
            dispatch(addTask(task));
        }
        const isDirectAccess = dataset.access === 'direct' || process.env.REACT_APP_STATIC === 'true';
        if (isDirectAccess) {
            dataset.api = new DirectAccessDataset();
        } else {
            dataset.api = new RestDataset();
        }

        const initPromise = dataset.api.init(id, dataset.url);
        const promises = [initPromise];

        promises.push(dataset.api.getJobs(id).then(jobs => {
            jobResults = jobs;
        }));

        const schemaPromise = dataset.api.getSchemaPromise().then(result => {
            newDataset = result;
        });
        promises.push(schemaPromise);
        if (!isDirectAccess) {
            const categoriesRenamePromise = getState().serverInfo.api.getCategoryNamesPromise(dataset.id).then(results => {
                categoryNameResults = results;
            });

            const savedViewsPromise = getState().serverInfo.api.getViewsPromise(dataset.id).then(results => {
                datasetViews = results;
            });
            promises.push(categoriesRenamePromise);
            promises.push(savedViewsPromise);
        }

        return Promise.all(promises).then(() => onPromisesComplete()).finally(() => {
            if (task) {
                dispatch(removeTask(task));
            }
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve dataset. Please try again.');
        });

    };
}


function handleSelectionResult(selectionResult, clear) {
    // not need to clear when adding a new feature
    return function (dispatch, getState) {
        const state = getState();
        if (selectionResult) {
            dispatch(setSelection(selectionResult.indices));
            // userPoints are in chart space, points are in server space, count is total number of cells selected
            if (selectionResult.summary) {
                // merge or clear selection
                let selectionSummary = clear ? selectionResult.summary : Object.assign({}, getState().featureSummary, selectionResult.summary);
                dispatch(setFeatureSummary(selectionSummary));
            }
            if (selectionResult.distribution) {
                let selectedDistributionData = state.selectedDistributionData;
                if (clear) {
                    selectedDistributionData = [];
                }
                const groupedSearchTokens = groupBy(state.searchTokens, 'type');
                addFeatureSetsToX(getFeatureSets(state.markers, groupedSearchTokens[FEATURE_TYPE.FEATURE_SET] || []), (groupedSearchTokens[FEATURE_TYPE.X] || []).map(item => item.id));
                selectedDistributionData = updateDistributionData(selectionResult.distribution, selectedDistributionData, groupedSearchTokens);
                dispatch(setSelectedDistributionData(selectedDistributionData));
            }
        }
    };
}

function updateDistributionData(newDistributionData, existingDistributionData, groupedSearchTokens) {
    let keys = new Set(Object.keys(existingDistributionData));
    if (newDistributionData) {
        for (const key in newDistributionData) {
            keys.add(key);
        }
    }
    keys.forEach(key => {
        existingDistributionData[key] = _updateDistributionData(newDistributionData ? newDistributionData[key] : null, existingDistributionData[key] || [], groupedSearchTokens);
    });
    return existingDistributionData;
}

function _updateDistributionData(newDistributionData, existingDistributionData, groupedSearchTokens) {
    const xTokens = (groupedSearchTokens[FEATURE_TYPE.X] || []).map(item => item.id);
    const obsCatTokens = (groupedSearchTokens[FEATURE_TYPE.OBS_CAT] || []).map(item => item.id);
    const obsTokens = (groupedSearchTokens[FEATURE_TYPE.OBS] || []).map(item => item.id);
    const moduleTokens = (groupedSearchTokens[FEATURE_TYPE.MODULE] || []).map(item => item.id);
    const dimensionKeys = [obsCatTokens.join('-')];
    // keep active dimensions and features only
    let distributionData = existingDistributionData.filter(entry => dimensionKeys.indexOf(entry.dimension) !== -1 &&
        (xTokens.indexOf(entry.feature) !== -1 || obsTokens.indexOf(entry.feature) !== -1 || moduleTokens.indexOf(entry.feature) !== -1));

    if (newDistributionData) {
        // remove old data that is also in new data
        newDistributionData.forEach(entry => {
            for (let i = 0; i < distributionData.length; i++) {
                if (distributionData[i].dimension === entry.dimension && distributionData[i].feature === entry.feature) {
                    distributionData.splice(i, 1);
                    break;
                }
            }
        });
        distributionData = distributionData.concat(newDistributionData);
    }

    // sort features matching X entry order
    const featureSortOrder = {};
    xTokens.concat(obsTokens).concat(moduleTokens).forEach((name, index) => {
        featureSortOrder[name] = index;
    });
    distributionData.sort((a, b) => {
        a = featureSortOrder[a.feature];
        b = featureSortOrder[b.feature];
        return a - b;
    });
    return distributionData;
}


function _updateCharts(onError, updateActiveFeature = true) {
    return function (dispatch, getState) {

        const state = getState();
        if (state.dataset == null) {
            return;
        }

        const groupedSearchTokens = groupBy(state.searchTokens, 'type');
        Object.values(FEATURE_TYPE).forEach(key => {
            if (!groupedSearchTokens[key]) {
                groupedSearchTokens[key] = []; // default
            }
        });
        const featureSetTokens = groupedSearchTokens[FEATURE_TYPE.FEATURE_SET];
        const xValues = groupedSearchTokens[FEATURE_TYPE.X].map(token => token.id);
        const obsCatValues = groupedSearchTokens[FEATURE_TYPE.OBS_CAT].map(token => token.id);
        const obsValues = groupedSearchTokens[FEATURE_TYPE.OBS].map(token => token.id);
        const moduleValues = groupedSearchTokens[FEATURE_TYPE.MODULE].map(token => token.id);

        addFeatureSetsToX(getFeatureSets(state.markers, featureSetTokens), xValues);

        const embeddings = state.embeddings;
        let distributionData = state.distributionData;
        let selectedDistributionData = state.selectedDistributionData;
        let embeddingData = state.embeddingData;
        const globalFeatureSummary = state.globalFeatureSummary;
        const cachedData = state.cachedData;
        const embeddingsToFetch = [];
        const valuesToFetch = new Set();
        const embeddingKeys = new Set();
        const features = new Set();
        const filterJson = getFilterJson(state);
        xValues.concat(obsValues).concat(obsCatValues).concat(moduleValues).forEach(feature => {
            features.add(feature);
        });
        const embeddingImagesToFetch = [];
        embeddings.forEach(embedding => {
            const embeddingKey = getEmbeddingKey(embedding);
            embeddingKeys.add(embeddingKey);
            if (cachedData[embeddingKey] == null) {
                if (embedding.dimensions > 0) {
                    embeddingsToFetch.push(embedding);
                } else {
                    if (cachedData[embedding.attrs.group] == null) {
                        valuesToFetch.add(embedding.attrs.group);
                    }

                    embedding.attrs.selection.forEach(selection => {
                        if (cachedData[selection[0]] == null) {
                            valuesToFetch.add(selection[0]);
                        }
                    });
                    embeddingImagesToFetch.push(embedding);
                }
            }
        });
        const promises = [];
        embeddingImagesToFetch.forEach(embedding => {
            const embeddingKey = getEmbeddingKey(embedding);
            if (cachedData[embeddingKey] == null) {
                const url = state.dataset.api.getFileUrl(embedding.image);
                const imagePromise = new Promise((resolve, reject) => {
                    fetch(url).then(result => result.text()).then(text => document.createRange().createContextualFragment(text).firstElementChild).then(node => {
                        // inline css
                        if (node.querySelector('style')) {
                            const div = document.createElement('div');
                            div.style.display = 'none';
                            div.appendChild(node);
                            document.body.appendChild(div);
                            const style = node.querySelector('style');
                            const rules = style.sheet.rules;

                            for (let i = 0; i < rules.length; i++) {
                                const rule = rules[i];
                                const matches = node.querySelectorAll(rule.selectorText);
                                const styleMap = rule.styleMap;
                                for (let j = 0; j < matches.length; j++) {
                                    const child = matches[j];
                                    for (let key of styleMap.keys()) {
                                        child.style[key] = styleMap.get(key).toString();
                                    }
                                }
                            }
                            div.remove();
                            style.remove();
                        }
                        cachedData[embeddingKey] = node;
                        resolve();
                    });
                });
                promises.push(imagePromise);
            }

        });
        if (filterJson != null) {
            let filterDependencies = getDatasetFilterDependencies(state.datasetFilter);
            filterDependencies.features.forEach(feature => {
                if (cachedData[feature] == null) {
                    valuesToFetch.add(feature);
                }
            });
        }
        features.forEach(feature => {
            if (cachedData[feature] == null) {
                valuesToFetch.add(feature);
            }
        });
        // don't fetch "other" features but add to features so they are displayed in embeddings
        // "other" features need to be in cachedData and globalFeatureSummary
        const keys = Object.values(FEATURE_TYPE);
        const otherSearchTokenKeys = [];
        for (const key in groupedSearchTokens) {
            if (keys.indexOf(key) === -1) {
                otherSearchTokenKeys.push(key);
            }
        }

        otherSearchTokenKeys.forEach(key => {
            groupedSearchTokens[key].forEach(item => features.add(item.id));
        });
        // set active flag on cached embedding data
        embeddingData.forEach(traceInfo => {
            const embeddingKey = getEmbeddingKey(traceInfo.embedding);
            const active = embeddingKeys.has(embeddingKey) && (features.has(traceInfo.name) || (features.size === 0 && traceInfo.name === '__count'));
            if (active) {
                traceInfo.date = new Date();
            }
            traceInfo.active = active;
        });


        const distributionCategories = obsCatValues.slice();
        const distributionCategoryKeys = [distributionCategories.join('-')];
        const distributionMeasuresToFetch = new Set();
        const distribution = (xValues.length > 0 || moduleValues.length > 0 || obsValues.length > 0 || otherSearchTokenKeys.length > 0) && obsCatValues.length > 0;
        if (distribution) { // TODO cleanup this code
            let cachedDistributionKeys = {}; // category-feature
            for (const key in distributionData) {
                distributionData[key].forEach(distributionDataItem => {
                    cachedDistributionKeys[distributionDataItem.name + '-' + distributionDataItem.feature] = true;
                });
            }
            distributionCategoryKeys.forEach(category => {
                xValues.forEach(feature => {
                    let key = category + '-' + feature;
                    if (cachedDistributionKeys[key] == null) {
                        distributionMeasuresToFetch.add(feature);
                    }
                });
                obsValues.forEach(feature => {
                    let key = category + '-' + feature;
                    if (cachedDistributionKeys[key] == null) {
                        distributionMeasuresToFetch.add('obs/' + feature);
                    }
                });
                moduleValues.forEach(feature => {
                    let key = category + '-' + feature;
                    if (cachedDistributionKeys[key] == null) {
                        distributionMeasuresToFetch.add('modules/' + feature);
                    }
                });
                otherSearchTokenKeys.forEach(searchTokenKey => {
                    groupedSearchTokens[searchTokenKey].forEach(item => {
                        let key = category + '-' + item.id;
                        if (cachedDistributionKeys[key] == null) {
                            distributionMeasuresToFetch.add(searchTokenKey + '/' + item.id);
                        }
                    });
                });
            });

        }
        let q = {};
        if (embeddingsToFetch.length > 0) {
            q.embedding = [];
            embeddingsToFetch.forEach(embedding => {
                if (embedding.mode != null) {// fetch unbinned coordinates
                    const key = getEmbeddingKey(embedding, false);
                    if (indexOf(embeddingsToFetch, (e => getEmbeddingKey(e) === key)) !== -1) {
                        const index = indexOf(state.dataset.embeddings, (e => getEmbeddingKey(e, false) === key));
                        if (index === -1) {
                            throw new Error(key + ' not found in ' + state.dataset.embeddings);
                        }
                        embeddingsToFetch.push(state.dataset.embeddings[index]);
                    }
                }
            });
            embeddingsToFetch.forEach(embedding => {
                q.embedding.push(embedding);
            });
        }
        if (valuesToFetch.size > 0) {
            const dataset = state.dataset;
            q.values = {measures: [], dimensions: []};
            valuesToFetch.forEach(value => {
                if (dataset.obsCat.indexOf(value) !== -1) {
                    q.values.dimensions.push(value);
                } else if (dataset.obs.indexOf(value) !== -1) {
                    q.values.measures.push('obs/' + value);
                } else if (moduleValues.indexOf(value) !== -1) {
                    q.values.measures.push('module/' + value);
                } else {
                    q.values.measures.push(value);
                }
            });
        }


        const globalFeatureSummaryMeasuresCacheMiss = [];
        const globalFeatureSummaryDimensionsCacheMiss = [];
        xValues.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryMeasuresCacheMiss.push(feature);
            }
        });
        obsValues.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryMeasuresCacheMiss.push('obs/' + feature);
            }
        });
        moduleValues.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryMeasuresCacheMiss.push('module/' + feature);
            }
        });
        obsCatValues.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryDimensionsCacheMiss.push(feature);
            }
        });


        if (globalFeatureSummaryMeasuresCacheMiss.length > 0 || globalFeatureSummaryDimensionsCacheMiss.length > 0) {
            q.stats = {
                measures: globalFeatureSummaryMeasuresCacheMiss,
                dimensions: globalFeatureSummaryDimensionsCacheMiss
            };
        }

        if (distributionCategories.length > 0 && distributionMeasuresToFetch.size > 0) {
            q.groupedStats = {
                measures: Array.from(distributionMeasuresToFetch),
                dimensions: [distributionCategories]
            };
        }


        // TODO update selection in new embedding space
        if (filterJson != null && (globalFeatureSummaryMeasuresCacheMiss.length > 0 || globalFeatureSummaryDimensionsCacheMiss.length > 0)) {
            q.selection = {
                filter: filterJson,
                measures: globalFeatureSummaryDimensionsCacheMiss.length > 0 ? xValues : globalFeatureSummaryMeasuresCacheMiss,
                dimensions: obsCatValues
            };
        }

        const fetchData = Object.keys(q).length > 0;
        const task = fetchData || promises.length > 0 ? {
            name: 'Update charts'
        } : null;
        if (task) {
            dispatch(addTask(task));
        }
        const dataPromise = fetchData ? state.dataset.api.getDataPromise(q, cachedData) : Promise.resolve({});
        const allPromises = [dataPromise].concat(promises);
        return Promise.all(allPromises).then(values => {
            const result = values[0];
            dispatch(setGlobalFeatureSummary(result.summary));
            const newEmbeddingData = getNewEmbeddingData(state, features);
            embeddingData = embeddingData.concat(newEmbeddingData);

            dispatch(setDistributionData(updateDistributionData(result.distribution, distributionData, groupedSearchTokens)));
            dispatch(setEmbeddingData(embeddingData));
            if (updateActiveFeature) {
                dispatch(setActiveFeature(getNewActiveFeature(embeddingData)));
            }
            if (result.selection) {
                dispatch(handleSelectionResult(result.selection, false));
            } else { // clear selection
                dispatch(setSelectedDistributionData(updateDistributionData(null, selectedDistributionData, groupedSearchTokens)));
            }
        }).finally(() => {
            if (task) {
                dispatch(removeTask(task));
            }
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve data. Please try again.');
            if (onError) {
                onError(err);
            }
        });
    };

}


function getNewActiveFeature(embeddingData) {
    let traces = embeddingData.filter(trace => trace.active);
    if (traces.length === 0) {
        return null;
    }
    const trace = traces[traces.length - 1];
    const embeddingKey = getTraceKey(trace); // last feature becomes primary
    return {name: trace.name, type: trace.featureType, embeddingKey: embeddingKey};
}


// depends on global feature summary
function getNewEmbeddingData(state, features) {
    const embeddings = state.embeddings;
    const embeddingData = state.embeddingData;
    const newEmbeddingData = [];
    const globalFeatureSummary = state.globalFeatureSummary;
    const interpolator = state.interpolator;
    const dataset = state.dataset;
    const cachedData = state.cachedData;
    const selection = state.selection;
    const searchTokens = state.searchTokens;
    const categoricalNames = state.categoricalNames;
    const existingFeaturePlusEmbeddingKeys = new Set();
    embeddingData.forEach(embeddingDatum => {
        const embeddingKey = getEmbeddingKey(embeddingDatum.embedding);
        const key = embeddingDatum.name + '_' + embeddingKey;
        existingFeaturePlusEmbeddingKeys.add(key);
    });
    if (features.size === 0) {
        features.add('__count');
    }

    embeddings.forEach(embedding => {
        const embeddingKey = getEmbeddingKey(embedding);
        // type can be image, scatter, or meta_image
        const traceType = embedding.spatial != null ? embedding.spatial.type : (embedding.type ? embedding.type : TRACE_TYPE_SCATTER);
        let coordinates = traceType !== TRACE_TYPE_META_IMAGE ? cachedData[embeddingKey] : null;
        if (coordinates == null && embedding.mode != null) {
            const unbinnedCoords = cachedData[getEmbeddingKey(embedding, false)];
            const binnedValues = createEmbeddingDensity(unbinnedCoords[embedding.name + '_1'], unbinnedCoords[embedding.name + '_2']);
            const binnedCoords = {};
            binnedCoords[embedding.name + '_1'] = binnedValues.x;
            binnedCoords[embedding.name + '_2'] = binnedValues.y;
            binnedCoords[embedding.name + '_index'] = binnedValues.index;
            binnedCoords[embedding.name + '_3'] = binnedValues.index;
            cachedData[embeddingKey] = binnedCoords; // save binned coords
            coordinates = binnedCoords;
        }
        const x = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_1'] : null;
        const y = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_2'] : null;
        const z = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_3'] : null;

        features.forEach(feature => {
            const featurePlusEmbeddingKey = feature + '_' + embeddingKey;
            let featureKey = feature;
            if (!existingFeaturePlusEmbeddingKeys.has(featurePlusEmbeddingKey)) {
                let featureSummary = globalFeatureSummary[feature];
                let values = cachedData[featureKey];
                if (values == null) {
                    if (feature === '__count') {
                        values = new Int8Array(dataset.shape[0]);
                        values.fill(1);
                    } else {
                        console.log(featureKey + ' not found');
                    }
                }
                let purity = null;
                if (values.value !== undefined) {
                    purity = values.purity;
                    values = values.value;
                }
                const searchToken = find(searchTokens, (item => item.id === feature));
                let featureType;
                if (searchToken) {
                    featureType = searchToken.type;
                }
                // could also be feature in a set
                if (featureType == null) {
                    featureType = feature === '__count' ? FEATURE_TYPE.COUNT : FEATURE_TYPE.X;
                }

                const isCategorical = featureType === FEATURE_TYPE.OBS_CAT;
                let colorScale = null;

                if (!isCategorical) {
                    // __count range is per embedding so we always recompute for now
                    if (isPlainObject(values)) {
                        let newValues = new Float32Array(dataset.shape[0]);
                        for (let i = 0, n = values.index.length; i < n; i++) {
                            newValues[values.index[i]] = values.values[i];
                        }
                        values = newValues;
                    }
                    if (featureSummary == null || feature === '__count') {

                        let min = Number.MAX_VALUE;
                        let max = -Number.MAX_VALUE;
                        let sum = 0;
                        for (let i = 0, n = values.length; i < n; i++) {
                            let value = values[i];
                            min = value < min ? value : min;
                            max = value > max ? value : max;
                            sum += value;
                        }
                        featureSummary = {min: min, max: max, mean: sum / values.length};
                        globalFeatureSummary[feature] = featureSummary;
                    }
                    let domain = [featureSummary.min, featureSummary.max];
                    if (featureSummary.customMin != null && !isNaN(featureSummary.customMin)) {
                        domain[0] = featureSummary.customMin;
                    }
                    if (featureSummary.customMax != null && !isNaN(featureSummary.customMax)) {
                        domain[1] = featureSummary.customMax;
                    }
                    let typeInterpolator = interpolator[featureType];
                    if (typeInterpolator == null) {
                        typeInterpolator = DEFAULT_INTERPOLATORS[FEATURE_TYPE.X];
                        interpolator[featureType] = typeInterpolator;
                    }
                    colorScale = createColorScale(typeInterpolator).domain(domain);

                } else {
                    let traceUniqueValues = featureSummary.categories;
                    // if (traceUniqueValues.length === 1 && traceUniqueValues[0] === true) {
                    //     traceUniqueValues = traceUniqueValues.concat([null]);
                    //     traceSummary.categories = traceUniqueValues;
                    //     traceSummary.counts = traceSummary.counts.concat([state.dataset.shape[0] - traceSummary.counts[0]]);
                    // }

                    let colorMap = dataset.colors ? dataset.colors[feature] : null;
                    let colors;
                    if (colorMap == null) {
                        if (traceUniqueValues.length <= 10) {
                            colors = schemeCategory10;
                        } else if (traceUniqueValues.length <= 12) {
                            colors = schemePaired;
                        } else if (traceUniqueValues.length <= 20) {
                            colors = CATEGORY_20B;
                        } else {
                            colors = CATEGORY_20B.concat(CATEGORY_20C);
                        }
                    } else {
                        if (isArray(colorMap)) {
                            colors = colorMap;
                        } else {
                            colors = [];
                            traceUniqueValues.forEach((value, index) => {
                                let color = colorMap[value];
                                if (color == null) {
                                    color = schemeCategory10[index % schemeCategory10.length];
                                }
                                colors.push(color);
                            });
                        }
                    }

                    // load saved colors from database
                    // category -> originalValue -> {newValue, positiveMarkers, negativeMarkers, color}
                    const originalValueToData = categoricalNames[feature];
                    if (originalValueToData) {
                        for (const originalValue in originalValueToData) {
                            for (const originalValue in originalValueToData) {
                                const value = originalValueToData[originalValue];
                                if (value.color != null) {
                                    const index = traceUniqueValues.indexOf(originalValue);
                                    if (index !== -1) {
                                        colors[index] = value.color;
                                    }
                                }
                            }
                        }
                    }
                    colorScale = scaleOrdinal(colors).domain(traceUniqueValues);
                    colorScale.summary = featureSummary;
                }

                if (traceType === TRACE_TYPE_META_IMAGE && embedding.categoryToIndices == null) {
                    const groupBy = cachedData[embedding.attrs.group];
                    const categoryToIndices = {};
                    const passingIndices = getPassingFilterIndices(cachedData, {filters: embedding.attrs.selection});
                    if (passingIndices.size === 0) {
                        throw new Error('No passing indices found');
                    }
                    for (let index of passingIndices) {
                        const category = groupBy[index];
                        let indices = categoryToIndices[category];
                        if (indices === undefined) {
                            indices = [];
                            categoryToIndices[category] = indices;
                        }
                        indices.push(index);
                    }
                    embedding.categoryToIndices = categoryToIndices;
                }

                let chartData = {
                    embedding: Object.assign({}, embedding),
                    name: feature,
                    featureType: featureType,
                    x: x,
                    y: y,
                    z: z != null ? z : undefined,
                    dimensions: z != null ? 3 : 2,
                    date: new Date(),
                    active: true,
                    colorScale: colorScale,
                    continuous: !isCategorical,
                    isCategorical: isCategorical,
                    values: values, // for color
                    type: traceType
                };
                if (chartData.mode != null) {
                    chartData.index = coordinates[embedding.name + '_index'];
                    chartData._values = chartData.values;
                    chartData.values = summarizeDensity(chartData.values, chartData.index, selection, chartData.continuous ? 'max' : 'mode');
                }
                if (traceType === TRACE_TYPE_SCATTER) {
                    chartData.positions = getPositions(chartData);
                }
                if (traceType === TRACE_TYPE_META_IMAGE) {
                    const svg = cachedData[getEmbeddingKey(embedding)];
                    chartData.source = svg.cloneNode(true);
                    chartData.zscore = true;
                    chartData.gallerySource = svg.cloneNode(true);
                    chartData.categoryToIndices = embedding.categoryToIndices;

                    if (chartData.continuous) {
                        // compute mean and standard deviation
                        colorScale.domain([-3, 3]);
                        let mean = 0;
                        let count = 0;
                        for (let category in embedding.categoryToIndices) {
                            const indices = embedding.categoryToIndices[category];
                            for (let i = 0, n = indices.length; i < n; i++) {
                                mean += chartData.values[indices[i]];
                                count++;
                            }
                        }
                        mean = mean / count;
                        let sum = 0;
                        for (let category in embedding.categoryToIndices) {
                            const indices = embedding.categoryToIndices[category];
                            for (let i = 0, n = indices.length; i < n; i++) {
                                let diff = chartData.values[indices[i]] - mean;
                                diff = diff * diff;
                                sum += diff;
                            }
                        }
                        const n = count - 1;
                        const variance = sum / n;
                        chartData.mean = mean;
                        chartData.stdev = Math.sqrt(variance);
                    }
                    chartData.fullCategoryToStats = createCategoryToStats(chartData, new Set());
                    chartData.categoryToStats = state.selection.size != null && state.selection.size === 0 ? chartData.fullCategoryToStats : createCategoryToStats(chartData, state.selection);
                }
                updateTraceColors(chartData);

                if (traceType === TRACE_TYPE_IMAGE) {
                    // TODO cache image
                    chartData.indices = !isCategorical ? indexSort(values, true) : randomSeq(values.length);
                    const url = dataset.api.getFileUrl(embedding.spatial.image);
                    chartData.tileSource = new OpenSeadragon.ImageTileSource({
                        url: url,
                        buildPyramid: true,
                        crossOriginPolicy: "Anonymous"
                    });
                }

                newEmbeddingData.push(chartData);
            }
        });
    });
    return newEmbeddingData;

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
        const task = {name: 'List Datasets'};
        dispatch(addTask(task));
        return getState().serverInfo.api.getDatasetsPromise()
            .then(choices => {
                dispatch(_setDatasetChoices(choices));
            })
            .finally(() => dispatch(removeTask(task)))
            .catch(err => {
                handleError(dispatch, err, 'Unable to retrieve datasets. Please try again.');
            });
    };
}


export function getDatasetStateJson(state) {
    const {
        chartOptions,
        combineDatasetFilters,
        activeFeature,
        dataset,
        embeddingLabels,
        distributionPlotOptions,
        distributionPlotInterpolator,
        embeddings,
        searchTokens,
        datasetFilter,
        interpolator,
        jobResultId,
        markerOpacity,
        pointSize,
        unselectedMarkerOpacity,
        distributionData
    } = state;

    let json = {
        dataset: dataset.id,
        embeddings: embeddings
    };
    if (jobResultId != null) {
        json.jobId = jobResultId;
    }

    const scatterPlot = chartOptions.scatterPlot;
    if (json.embeddings.length > 0 && scatterPlot != null) {
        json.camera = scatterPlot.getCameraDef();

    }
    if (activeFeature != null) {
        json.activeFeature = activeFeature;
    }
    let jsonChartOptions = {};

    const defaultChartOptions = {
        showFog: DEFAULT_SHOW_FOG, darkMode: DEFAULT_DARK_MODE,
        labelFontSize: DEFAULT_LABEL_FONT_SIZE,
        labelStrokeWidth: DEFAULT_LABEL_STROKE_WIDTH
    };

    for (let key in defaultChartOptions) {
        let value = chartOptions[key];
        if (value !== defaultChartOptions[key]) {
            jsonChartOptions[key] = value;
        }
    }
    if (pointSize !== DEFAULT_POINT_SIZE) {
        json.pointSize = pointSize;
    }


    if (Object.keys(jsonChartOptions).length > 0) {
        json.chartOptions = jsonChartOptions;
    }

    if (searchTokens.length > 0) {
        json.q = searchTokens;
    }

    let datasetFilterJson = {};
    for (let key in datasetFilter) {
        let filterObject = datasetFilter[key];
        if (Array.isArray(filterObject)) { // brush filter
            const array = [];
            filterObject.forEach(brush => {
                const brushJson = Object.assign({}, brush);
                brushJson.indices = Array.from(brush.indices);
                array.push(brushJson);
            });
            datasetFilterJson[key] = array;
        } else if (filterObject.operation === 'in') {
            datasetFilterJson[key] = {operation: filterObject.operation, value: filterObject.value};
        } else {
            if (filterObject.operation !== '' && !isNaN(filterObject.value) && filterObject.value != null) {
                datasetFilterJson[key] = {operation: filterObject.operation, value: filterObject.value};
            }
        }
    }
    if (Object.keys(datasetFilterJson).length > 0) {
        json.datasetFilter = datasetFilterJson;
    }
    if (combineDatasetFilters !== 'and') {
        json.combineDatasetFilters = combineDatasetFilters;
    }
    if (markerOpacity !== DEFAULT_MARKER_OPACITY) {
        json.markerOpacity = markerOpacity;
    }
    if (unselectedMarkerOpacity !== DEFAULT_UNSELECTED_MARKER_OPACITY) {
        json.unselectedMarkerOpacity = unselectedMarkerOpacity;
    }

    if (distributionData && distributionData.length > 0) {
        json.distributionPlotOptions = distributionPlotOptions;
        json.distributionPlotInterpolator = Object.assign({}, distributionPlotInterpolator, {value: null});
    }

    // TODO save custom color ranges per feature
    for (let key in interpolator) {
        const typedInterpolator = interpolator[key];
        const defaultInterpolator = DEFAULT_INTERPOLATORS[key];

        if (defaultInterpolator == null || typedInterpolator.name !== defaultInterpolator.name || typedInterpolator.reversed !== defaultInterpolator.reversed) {
            if (json.interpolator == null) {
                json.interpolator = {};
            }
            json.interpolator[key] = Object.assign({}, typedInterpolator, {value: undefined});
        }
    }

    if (embeddingLabels.length > 0) {
        json.embeddingLabels = embeddingLabels;
    }
    return json;
}

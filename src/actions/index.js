import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10, schemePaired} from 'd3-scale-chromatic';
import {saveAs} from 'file-saver';
import {find, findIndex, isArray, isEqual} from 'lodash';
import OpenSeadragon from 'openseadragon';
import isPlainObject from 'react-redux/lib/utils/isPlainObject';
import CustomError from '../CustomError';
import {getPassingFilterIndices} from '../dataset_filter';
import {DirectAccessDataset} from '../DirectAccessDataset';
import {updateJob} from '../JobResultsPanel';
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
    FEATURE_TYPE,
    getFeatureSets,
    getInterpolator,
    indexSort,
    randomSeq,
    splitSearchTokens,
    TRACE_TYPE_IMAGE,
    TRACE_TYPE_META_IMAGE,
    TRACE_TYPE_SCATTER,
    updateTraceColors
} from '../util';

export const API = process.env.REACT_APP_API_URL || 'api';

const authScopes = [
    'email',
    // 'profile',
    // 'https://www.googleapis.com/auth/userinfo.profile',
    // 'https://www.googleapis.com/auth/contacts.readonly',
    // 'https://www.googleapis.com/auth/devstorage.full_control',
];
export const SET_EMBEDDING_LABELS = 'SET_EMBEDDING_LABELS';
export const SET_DISTRIBUTION_PLOT_INTERPOLATOR = 'SET_DISTRIBUTION_PLOT_INTERPOLATOR';
export const SET_CHART_OPTIONS = 'SET_CHART_OPTIONS';
export const SET_COMBINE_DATASET_FILTERS = 'SET_COMBINE_DATASET_FILTERS';
export const SET_DATASET_FILTERS = 'SET_DATASET_FILTERS'; // saved dataset filters
export const SET_UNSELECTED_MARKER_SIZE = 'SET_UNSELECTED_MARKER_SIZE';

export const SET_ACTIVE_FEATURE = 'SET_ACTIVE_FEATURE';
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

export const SET_SELECTION = 'SET_SELECTION';
export const SET_FEATURE_SUMMARY = 'SET_FEATURE_SUMMARY';
export const SET_SAVED_DATASET_FILTER = 'SET_SAVED_DATASET_FILTER';
export const SET_SEARCH_TOKENS = 'SET_SEARCH_TOKENS';

export const SET_SELECTED_EMBEDDING = 'SET_SELECTED_EMBEDDING';
export const SET_MESSAGE = 'SET_MESSAGE';
export const SET_INTERPOLATOR = 'SET_INTERPOLATOR';
export const SET_POINT_SIZE = 'SET_POINT_SIZE';

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

export const SET_DISTRIBUTION_DATA = 'SET_DISTRIBUTION_DATA';
export const SET_SELECTED_DISTRIBUTION_DATA = 'SET_SELECTED_DISTRIBUTION_DATA';

export const SET_DISTRIBUTION_PLOT_OPTIONS = 'SET_DISTRIBUTION_PLOT_OPTIONS';
export const SET_EMBEDDING_DATA = 'SET_EMBEDDING_DATA';

export const SET_LOADING = 'SET_LOADING';
export const SET_TAB = 'SET_TAB';

export const SET_LOADING_APP = 'LOADING_APP';
export const SET_CACHED_DATA = 'SET_CACHED_DATA';

export const SET_JOB_RESULTS = 'SET_JOB_RESULTS';
export const SET_JOB_RESULT = 'SET_JOB_RESULT';

function isEmbeddingBinned(embedding) {
    return embedding.bin || embedding.precomputed;
}

export function getEmbeddingKey(embedding) {
    let fullName = embedding.name;
    if (embedding.dimensions) {
        fullName += '_' + embedding.dimensions;
    }
    if (embedding.bin || embedding.precomputed) {
        fullName += '_' + embedding.nbins + '_' + embedding.agg;
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
        if (process.env.REACT_APP_STATIC === 'true') {
            const serverInfo = {clientId: '', dynamic: false};
            serverInfo.api = new StaticServerApi();
            dispatch(setServerInfo(serverInfo));
            dispatch(_setLoadingApp({loading: false}));
            dispatch(listDatasets()).then(() => {
                dispatch(_loadSavedView());
            });
            return Promise.resolve();
        }
        return fetch(API + '/server').then(result => result.json()).then(serverInfo => {
            serverInfo.api = new RestServerApi();
            if (serverInfo.clientId == null) {
                serverInfo.clientId = '';
            }
            if (serverInfo.dynamic == null) {
                serverInfo.dynamic = true;
            }
            dispatch(setServerInfo(serverInfo));
            if (serverInfo.clientId === '' || serverInfo.clientId == null) { // no login required

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
                    datasetFilter[field] = {operation: item[1], value: item[2]};

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


export function submitJob(jobData) {
    return function (dispatch, getState) {
        let jobId;
        let timeout = 5 * 1000; // TODO
        function getJobStatus() {
            getState().dataset.api.getJob(jobId, false)
                .then(result => {
                    const jobResult = find(getState().jobResults, item => item.id === jobId);
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
            continue;
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
        dispatch(_setLoading(true));
        const state = getState();
        let filter = getFilterJson(state, true);

        state.dataset.api.getSelectedIdsPromise({
            filter: filter
        }, state.cachedData).then(result => {
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


export function setDistributionPlotOptions(payload) {
    return {type: SET_DISTRIBUTION_PLOT_OPTIONS, payload: payload};
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

export function setActiveFeature(payload) {
    return {type: SET_ACTIVE_FEATURE, payload: payload};
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
        // const isCurrentSelectionEmpty = state.selection.chart == null || Object.keys(state.selection.chart).length === 0;

        if (filter == null) {
            // if (!isCurrentSelectionEmpty) {
            //     dispatch(setSelection({chart: {}}));
            // }
            if (state.selection.size !== 0) {
                dispatch(setSelection(new Set()));
            }
            dispatch(setSelectedDistributionData([]));
            dispatch(setFeatureSummary({}));
            dispatch(_setLoading(false));
            return;
        }
        // let selectedEmbeddings = state.embeddings;
        // q.selection.embeddings = selectedEmbeddings.map(e => {
        //     return getEmbeddingJson(e);
        // });
        const cachedData = state.cachedData;
        getState().dataset.api.getDataPromise(q, cachedData).then(result => {
            dispatch(handleSelectionResult(result.selection, true));
        }).catch(err => {
            handleError(dispatch, err);
        }).finally(() => dispatch(_setLoading(false)));
    };
}


export function handleBrushFilterUpdated(payload) {
    return function (dispatch, getState) {
        const name = payload.name; // full basis name
        const value = payload.value;  // value has basis and indices
        const clear = payload.clear;
        const indices = payload.value != null ? payload.value.indices : null;

        let datasetFilter = getState().datasetFilter;
        // if (value && isEmbeddingBinned(value.basis)) {
        //     const bins = getState().cachedData[getEmbeddingKey(value.basis)].bins;
        //     value.points = convertIndicesToBins(value.points, bins);
        // }
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
        const p = getState().serverInfo.dynamic ? getState().serverInfo.api.setCategoryNamePromise({
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


export function restoreView(payload) {
    return {type: RESTORE_VIEW, payload: payload};
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


export function setServerInfo(payload) {
    return {type: SET_SERVER_INFO, payload: payload};
}

function loadDefaultDatasetView() {
    return function (dispatch, getState) {
        const dataset = getState().dataset;
        const embeddingNames = dataset.embeddings.map(e => e.name);
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
            const selectedEmbedding = dataset.embeddings[dataset.embeddings.map(e => e.name).indexOf(embeddingName)];
            let obsCat = null;
            let catPriorities = ['leiden', 'louvain', 'seurat_clusters', 'clusters'];
            for (let priorityIndex = 0; priorityIndex < catPriorities.length && obsCat == null; priorityIndex++) {
                for (let i = 0; i < dataset.obsCat.length; i++) {
                    if (dataset.obsCat[i].toLowerCase().indexOf(catPriorities[priorityIndex]) !== -1) {
                        obsCat = dataset.obsCat[i];
                        break;
                    }
                }
            }
            if (obsCat != null) {
                dispatch(setSearchTokensDirectly([{value: obsCat, type: FEATURE_TYPE.OBS_CAT}]));
            }
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
        if (savedView.colorScheme != null) {
            savedView.colorScheme.value = getInterpolator(savedView.colorScheme.name);
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
                if (savedView.distributionPlotOptions != null) {
                    dispatch(setDistributionPlotOptions(savedView.distributionPlotOptions));
                }
                let activeFeature = savedView.activeFeature;
                if (activeFeature == null && savedView.embeddings && savedView.embeddings.length > 0 && savedView.q != null && savedView.q.length > 0) {
                    activeFeature = {
                        name: savedView.q[0],
                        type: getFeatureType(getState().dataset, savedView.q[0]),
                        embeddingKey: savedView.q[0] + '_' + getEmbeddingKey(savedView.embeddings[0])
                    };
                }
                if (activeFeature != null) {
                    dispatch(setActiveFeature(activeFeature));
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
        const file = payload.file;
        const readers = payload.readers;
        const description = payload.description;
        const title = payload.title;
        const species = payload.species;
        let existingDataset = payload.dataset;
        const isEdit = existingDataset != null;
        if (existingDataset == null) {
            existingDataset = {};
        }
        const data = {};
        if (existingDataset.readers == null) {
            data.readers = readers;
        } else {
            const existingReaders = new Set(existingDataset.readers);
            let setsEqual = existingReaders.size === readers.length;
            if (setsEqual) {
                for (let i = 0; i < readers.length; i++) {
                    if (!existingReaders.has(readers[i])) {
                        setsEqual = false;
                        break;
                    }
                }
            }
            if (!setsEqual) {
                data.readers = readers;
            }
        }
        if (name != null && name !== existingDataset.name) {
            data.name = name;
        }
        if (description != null && description !== existingDataset.description) {
            data.description = description;
        }
        if (title != null && title !== existingDataset.title) {
            data.title = title;
        }
        if (species != null && species !== existingDataset.species) {
            data.species = species;
        }
        if (url != null && url !== existingDataset.url) {
            data.url = url;
        }
        if (file != null) {
            data.file = file;
        }
        if (Object.keys(data).length === 0) {
            return;
        }
        if (isEdit) {
            data.id = payload.dataset.id;
        }
        dispatch(_setLoading(true));
        const request = getState().serverInfo.api.upsertDatasetPromise(data);
        request.upload.addEventListener('progress', function (e) {
            if (file) {
                let percent = (e.loaded / e.total) * 100;
                dispatch(setMessage('Percent ' + percent));
            }
        });
        request.addEventListener('load', function (e) {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
            const status = request.status;
            if (status != 200) {
                return dispatch(setMessage('Unable to save dataset. Please try again.'));
            }
            // request.response holds response from the server
            const resp = JSON.parse(request.response);
            data.id = resp.id;
            data.owner = true;

            if (isEdit) {
                dispatch(updateDataset(data));
                dispatch(setMessage('Dataset updated'));
            } else {
                dispatch(_addDataset(data));
                dispatch(setMessage('Dataset added'));
            }
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

export function setJobResults(payload) {
    return {type: SET_JOB_RESULTS, payload: payload};
}

function _setJobResult(payload) {
    return {type: SET_JOB_RESULT, payload: payload};
}

export function deleteJobResult(payload) {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        getState().dataset.api.deleteJob(payload).then(() => {
            let jobResults = getState().jobResults;
            const index = findIndex(jobResults, item => item.id === payload);
            jobResults.splice(index, 1);
            if (getState().jobResult === payload) {
                dispatch(_setJobResult(null));
            }
            dispatch(setJobResults(jobResults.slice()));
        }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete result. Please try again.');
        });
    };
}

export function setJobResult(payload) {
    return function (dispatch, getState) {
        let jobResults = getState().jobResults;
        let jobResult = find(jobResults, item => item.id === payload);
        if (jobResult.data != null) { // data already loaded
            updateJob(jobResult);
            return dispatch(_setJobResult(payload));
        }
        dispatch(_setLoading(true));
        getState().dataset.api.getJob(payload, true).then((result) => {
            let jobResults = getState().jobResults;
            const jobResult = find(jobResults, item => item.id === payload);
            for (let key in result) {
                jobResult[key] = result[key];
            }
            updateJob(jobResult);
            dispatch(_setJobResult(payload));
        }).finally(() => {
            dispatch(_setLoading(false));
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
        let embeddingLabels = getState().embeddingLabels;
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


/**
 *
 * @param value
 * @param type one of X, obs, featureSet
 */
export function setSearchTokens(values, type) {
    return function (dispatch, getState) {
        const state = getState();
        let searchTokens = state.searchTokens;
        // keep all other types
        let removeType = [type];
        if (type === FEATURE_TYPE.OBS) {
            removeType.push(FEATURE_TYPE.OBS_CAT);
        }
        searchTokens = searchTokens.filter(item => removeType.indexOf(item.type) === -1);
        if (type === FEATURE_TYPE.OBS) {
            const obsCat = state.dataset.obsCat;
            values.forEach(val => {
                const type = obsCat.indexOf(val) !== -1 ? FEATURE_TYPE.OBS_CAT : FEATURE_TYPE.OBS;
                searchTokens.push({value: val, type: type});
            });
        } else {
            searchTokens = searchTokens.concat(values.map(item => {
                return {value: item, type: type};
            }));
        }
        dispatch({type: SET_SEARCH_TOKENS, payload: searchTokens.slice()});
        dispatch(_updateCharts());
    };
}

export function setSearchTokensDirectly(tokens) {
    return function (dispatch, getState) {
        dispatch({type: SET_SEARCH_TOKENS, payload: tokens});
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

        let categoryNameResults;
        let datasetFilters = [];
        let newDataset;
        let jobResults = [];

        function onPromisesComplete() {
            newDataset = Object.assign({}, dataset, newDataset);
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

            if (newDataset.results && newDataset.results.length > 0) {
                jobResults = jobResults.concat(newDataset.results);
            }
            dispatch(setJobResults(jobResults));
            // if (jobResults.length === 1) {
            //     dispatch(setJobResult(jobResults[0].id));
            // }

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
                dispatch(loadDefaultDatasetView());
            }
        }

        const isDirectAccess = dataset.access === 'direct';
        if (isDirectAccess) {
            dataset.api = new DirectAccessDataset();
        } else {
            dataset.api = new RestDataset();
        }

        const initPromise = dataset.api.init(id, dataset.url);
        const jobsPromise = dataset.api.getJobs(id).then(jobs => {
            jobResults = jobs;
        });
        const schemaPromise = dataset.api.getSchemaPromise().then(result => {
            newDataset = result;
        });
        const promises = [initPromise, jobsPromise, schemaPromise];
        if (!isDirectAccess) {
            const categoriesRenamePromise = getState().serverInfo.api.getCategoryNamesPromise(dataset.id).then(results => {
                categoryNameResults = results;
            });

            const filtersPromise = getState().serverInfo.api.getFiltersPromise(dataset.id).then(results => {
                datasetFilters = results;
            });
            promises.push(categoriesRenamePromise);
            promises.push(filtersPromise);
        }

        return Promise.all(promises).then(() => onPromisesComplete()).finally(() => {
            if (setLoading) {
                dispatch(_setLoading(false));
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
            // const indices = selectionResult.indices;
            // if (indices != null) {
            // let chartSelection = {};
            // for (let key in coordinates) {
            //     const selectedIndicesOrBins = coordinates[key].indices_or_bins;
            //     const embedding = state.embeddings[state.embeddings.map(e => getEmbeddingKey(e)).indexOf(key)];
            //     if (embedding == null) {
            //         console.log(key + ' missing from coordinates');
            //         continue;
            //     }
            //     let selectedPoints = selectedIndicesOrBins;
            //     if (isEmbeddingBinned(embedding)) {
            //         const coords = state.cachedData[getEmbeddingKey(embedding)];
            //         const bins = coords.bins;
            //         selectedPoints = convertBinsToIndices(bins, selectedIndicesOrBins);
            //     }
            //
            //     chartSelection[key] = {
            //         userPoints: new Set(selectedPoints),
            //         points: selectedIndicesOrBins
            //     };
            // }

            // } else {
            //     const isCurrentSelectionEmpty = state.selection.indices == null || Object.keys(state.selection.chart).length === 0;
            //     if (clear && !isCurrentSelectionEmpty) {
            //         dispatch(setSelection({count: selectionResult.count, indices: null}));
            //     }
            // }

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
                const searchTokens = splitSearchTokens(state.searchTokens);
                addFeatureSetsToX(getFeatureSets(state.markers, searchTokens.featureSets), searchTokens.X);
                selectedDistributionData = updateDistributionData(selectionResult.distribution, selectedDistributionData, searchTokens);
                dispatch(setSelectedDistributionData(selectedDistributionData));
            }
        }
    };
}

function updateDistributionData(newDistributionData, distributionData, searchTokens) {
    let dimensionKeys = [searchTokens.obsCat.join('-')];
    // keep active dimensions and features only
    distributionData = distributionData.filter(entry => dimensionKeys.indexOf(entry.dimension) !== -1 && searchTokens.X.indexOf(entry.feature) !== -1);

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
    let featureSortOrder = {};
    searchTokens.X.forEach((name, index) => {
        featureSortOrder[name] = index;
    });
    distributionData.sort((a, b) => {
        a = featureSortOrder[a.feature];
        b = featureSortOrder[b.feature];
        return a - b;
    });
    return distributionData;
}

function _updateCharts(onError) {
    return function (dispatch, getState) {
        const state = getState();
        if (state.dataset == null) {
            return;
        }
        dispatch(_setLoading(true));
        const searchTokens = splitSearchTokens(state.searchTokens);
        addFeatureSetsToX(getFeatureSets(state.markers, searchTokens.featureSets), searchTokens.X);
        let distribution = (searchTokens.X.length > 0 || searchTokens.featureSets.length > 0) && searchTokens.obsCat.length > 0;
        const embeddings = state.embeddings;
        let distributionData = state.distributionData;
        let selectedDistributionData = state.selectedDistributionData;
        let embeddingData = state.embeddingData;
        const globalFeatureSummary = state.globalFeatureSummary;
        const cachedData = state.cachedData;
        const embeddingCoordinatesToFetch = [];
        const valuesToFetch = new Set();
        const binnedEmbeddingValuesToFetch = [];
        const embeddingKeys = new Set();
        const features = new Set();
        const filterJson = getFilterJson(state);
        searchTokens.X.concat(searchTokens.obs).concat(searchTokens.obsCat).concat(searchTokens.metafeatures).forEach(feature => {
            features.add(feature);
        });
        const embeddingImagesToFetch = [];
        embeddings.forEach(selectedEmbedding => {
            const embeddingKey = getEmbeddingKey(selectedEmbedding);
            embeddingKeys.add(embeddingKey);

            if (cachedData[embeddingKey] == null) {
                if (selectedEmbedding.dimensions > 0) {
                    embeddingCoordinatesToFetch.push(selectedEmbedding);
                } else {
                    if (cachedData[selectedEmbedding.attrs.group] == null) {
                        valuesToFetch.add(selectedEmbedding.attrs.group);
                    }

                    selectedEmbedding.attrs.selection.forEach(selection => {
                        if (cachedData[selection[0]] == null) {
                            valuesToFetch.add(selection[0]);
                        }
                    });
                    embeddingImagesToFetch.push(selectedEmbedding);
                }
            }

            if (isEmbeddingBinned(selectedEmbedding)) {
                const data = {values: [], embedding: selectedEmbedding};
                features.forEach(feature => {
                    let key = feature + '_' + embeddingKey;
                    if (cachedData[key] == null) {
                        data.values.push(feature);
                    }
                });
                if (data.values.length > 0) {
                    binnedEmbeddingValuesToFetch.push(data);
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
        // set active flag on cached embedding data
        embeddingData.forEach(traceInfo => {
            const embeddingKey = getEmbeddingKey(traceInfo.embedding);
            const active = embeddingKeys.has(embeddingKey) && (features.has(traceInfo.name) || (features.size === 0 && traceInfo.name === '__count'));
            if (active) {
                traceInfo.date = new Date();
            }
            traceInfo.active = active;
        });

        let distributionCategories = searchTokens.obsCat.slice();
        let distributionCategoryKeys = [distributionCategories.join('-')];
        let distributionMeasuresToFetch = new Set();
        if (distribution) {
            let cachedDistributionKeys = {}; // category-feature
            distributionData.forEach(distributionDataItem => {
                cachedDistributionKeys[distributionDataItem.name + '-' + distributionDataItem.feature] = true;
            });
            distributionCategoryKeys.forEach(category => {
                searchTokens.X.forEach(feature => {
                    let key = category + '-' + feature;
                    if (cachedDistributionKeys[key] == null) {
                        distributionMeasuresToFetch.add(feature);
                    }
                });
            });
        }
        let q = {};
        if (embeddingCoordinatesToFetch.length > 0 || binnedEmbeddingValuesToFetch.length > 0) {
            q.embedding = [];

            const embeddingCoordinatesToFetchKeys = embeddingCoordinatesToFetch.map(e => getEmbeddingKey(e));
            binnedEmbeddingValuesToFetch.forEach(datum => {
                const embeddingJson = getEmbeddingJson(datum.embedding);
                embeddingJson.dimensions = [];
                embeddingJson.measures = [];
                embeddingJson.coords = embeddingCoordinatesToFetchKeys.indexOf(getEmbeddingKey(datum.embedding)) !== -1;
                q.embedding.push(embeddingJson);
                datum.values.forEach(value => {
                    if (searchTokens.obsCat.indexOf(value) !== -1) {
                        embeddingJson.dimensions.push(value);
                    } else if (searchTokens.obs.indexOf(value) !== -1) {
                        embeddingJson.measures.push('obs/' + value);
                    } else {
                        embeddingJson.measures.push(value);
                    }
                });
            });
            const binnedEmbeddingValuesToFetchKeys = binnedEmbeddingValuesToFetch.map(e => getEmbeddingKey(e.embedding));
            embeddingCoordinatesToFetch.forEach(embedding => {
                if (binnedEmbeddingValuesToFetchKeys.indexOf(getEmbeddingKey(embedding)) === -1) {
                    const embeddingJson = getEmbeddingJson(embedding);
                    embeddingJson.coords = true;
                    q.embedding.push(embeddingJson);
                }
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
                } else if (searchTokens.metafeatures.indexOf(value) !== -1) {
                    q.values.measures.push('metagenes/' + value);
                } else {
                    q.values.measures.push(value);
                }
            });
        }

        let globalFeatureSummaryXCacheMiss = [];
        let globalFeatureSummaryObsContinuousCacheMiss = [];
        let globalFeatureSummaryDimensionsCacheMiss = [];
        searchTokens.X.forEach(feature => {
            if (globalFeatureSummary[feature] == null) {
                globalFeatureSummaryXCacheMiss.push(feature);
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

        if (distributionCategories.length > 0 && distributionMeasuresToFetch.size > 0) {
            q.groupedStats = {
                measures: Array.from(distributionMeasuresToFetch),
                dimensions: [distributionCategories]
            };
        }


        // TODO update selection in new embedding space
        if (filterJson != null && (globalFeatureSummaryXCacheMiss.length > 0 || globalFeatureSummaryDimensionsCacheMiss.length > 0)) {
            q.selection = {
                filter: filterJson,
                measures: globalFeatureSummaryDimensionsCacheMiss.length > 0 ? searchTokens.X : globalFeatureSummaryXCacheMiss,
                dimensions: searchTokens.obsCat
            };
        }

        let dataPromise = Object.keys(q).length > 0 ? state.dataset.api.getDataPromise(q, cachedData) : Promise.resolve({});
        const allPromises = [dataPromise].concat(promises);
        return Promise.all(allPromises).then(values => {
            const result = values[0];
            const newSummary = result.summary || {};
            for (let key in newSummary) {
                globalFeatureSummary[key] = newSummary[key];
            }
            dispatch(setGlobalFeatureSummary(globalFeatureSummary));
            updateEmbeddingData(state, features);
            dispatch(setDistributionData(updateDistributionData(result.distribution, distributionData, searchTokens)));
            dispatch(setEmbeddingData(embeddingData.slice()));
            if (result.selection) {
                dispatch(handleSelectionResult(result.selection, false));
            } else { // clear selection
                dispatch(setSelectedDistributionData(updateDistributionData(null, selectedDistributionData, searchTokens)));
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


function getFeatureType(dataset, feature) {
    if (feature === '__count') {
        return FEATURE_TYPE.COUNT;
    } else if (dataset.obsCat.indexOf(feature) !== -1) {
        return FEATURE_TYPE.OBS_CAT;
    } else if (dataset.obs.indexOf(feature) !== -1) {
        return FEATURE_TYPE.OBS;
    }
    return FEATURE_TYPE.X;
}

// depends on global feature summary
function updateEmbeddingData(state, features) {
    const embeddings = state.embeddings;
    let embeddingData = state.embeddingData;
    const globalFeatureSummary = state.globalFeatureSummary;
    const interpolator = state.interpolator;
    const dataset = state.dataset;
    const cachedData = state.cachedData;
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
        const isBinned = isEmbeddingBinned(embedding);
        // type can be image, scatter, or meta_image
        const traceType = embedding.spatial != null ? embedding.spatial.type : (embedding.type ? embedding.type : TRACE_TYPE_SCATTER);
        const coordinates = traceType !== TRACE_TYPE_META_IMAGE ? cachedData[embeddingKey] : null;
        const x = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_1'] : null;
        const y = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_2'] : null;
        const z = traceType !== TRACE_TYPE_META_IMAGE ? coordinates[embedding.name + '_3'] : null;


        const bins = isBinned ? coordinates.bins : null;
        features.forEach(feature => {
            const featurePlusEmbeddingKey = feature + '_' + embeddingKey;
            let featureKey = feature;
            if (!existingFeaturePlusEmbeddingKeys.has(featurePlusEmbeddingKey)) {

                if (isBinned) {
                    featureKey = featurePlusEmbeddingKey;
                }
                let featureSummary = globalFeatureSummary[feature];
                let values = cachedData[featureKey];
                if (values == null && feature === '__count' && !isBinned) {
                    values = new Int8Array(dataset.shape[0]);
                    values.fill(1);
                }
                let purity = null;

                if (values.value !== undefined) {
                    purity = values.purity;
                    values = values.value;
                }
                const featureType = getFeatureType(dataset, feature);
                let isCategorical = featureType === FEATURE_TYPE.OBS_CAT;
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
                    colorScale = createColorScale(interpolator).domain(domain);

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
                    bins: bins,
                    dimensions: z != null ? 3 : 2,
                    npoints: values.length,
                    date: new Date(),
                    active: true,
                    colorScale: colorScale,
                    continuous: !isCategorical,
                    isCategorical: isCategorical,
                    values: values, // for color
                    type: traceType
                    // purity: purity,
                    // text: values,
                };
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
                    chartData.categoryToStats = state.selection.size === 0 ? chartData.fullCategoryToStats : createCategoryToStats(chartData, state.selection);
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

                embeddingData.push(chartData);
            }
        });
    });

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



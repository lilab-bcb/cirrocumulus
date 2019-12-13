import {color} from 'd3-color';
import {scaleLinear, scaleOrdinal, scaleSequential} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import {saveAs} from 'file-saver';
import {cloneDeep, isEqual, isPlainObject} from 'lodash';
import CustomError from '../CustomError';
import PlotUtil, {CATEGORY_20B, CATEGORY_20C, getInterpolator} from '../PlotUtil';

//export const API = 'http://localhost:5000/api';
export const API = '/api';

const authScopes = [
    'email',
    'profile',
    'https://www.googleapis.com/auth/userinfo.profile',
    'https://www.googleapis.com/auth/contacts.readonly',
    'https://www.googleapis.com/auth/devstorage.full_control',
];

export const SET_DATASET_FILTERS = 'SET_DATASET_FILTERS'; // saved dataset filters
export const SET_UNSELECTED_MARKER_SIZE = 'SET_UNSELECTED_MARKER_SIZE';
export const SET_UNSELECTED_MARKER_SIZE_UI = 'SET_UNSELECTED_MARKER_SIZE_UI';
export const SET_SERVER_INFO = "SET_SERVER_INFO";
export const SET_DATASET_FILTER = 'SET_DATASET_FILTER';
export const ADD_DATASET = 'ADD_DATASET';
export const DELETE_DATASET = 'DELETE_DATASET';
export const UPDATE_DATASET = 'UPDATE_DATASET';
export const SET_GLOBAL_FEATURE_SUMMARY = 'SET_GLOBAL_FEATURE_SUMMARY';

export const SET_MARKER_SIZE = 'SET_MARKER_SIZE';
export const SET_MARKER_OPACITY = 'SET_MARKER_OPACITY';

export const SET_EMBEDDING_CHART_SIZE = "SET_EMBEDDING_CHART_SIZE";

export const SET_UNSELECTED_MARKER_OPACITY = 'SET_UNSELECTED_MARKER_OPACITY';

// update chart

export const SET_SELECTION = 'SET_SELECTION';
export const SET_FEATURE_SUMMARY = 'SET_FEATURE_SUMMARY';
export const SET_SAVED_DATASET_FILTER = 'SET_SAVED_DATASET_FILTER';
export const SET_FEATURES = 'SET_FEATURES';
export const SET_FEATURES_UI = 'SET_FEATURES_UI';
export const SET_GROUP_BY = 'SET_GROUP_BY';
export const SET_SELECTED_EMBEDDING = 'SET_SELECTED_EMBEDDING';
export const SET_NUMBER_OF_BINS = 'SET_NUMBER_OF_BINS';
export const SET_BIN_VALUES = 'SET_BIN_VALUES';
export const SET_BIN_SUMMARY = 'SET_BIN_SUMMARY';
export const SET_MESSAGE = 'SET_MESSAGE';
export const SET_INTERPOLATOR = 'SET_INTERPOLATOR';

// update ui
export const SET_EMAIL = 'SET_EMAIL';
export const SET_USER = 'SET_USER';
export const SET_DATASET = 'SET_DATASET';

export const SET_DIALOG = 'SET_DIALOG';

export const EDIT_DATASET_DIALOG = 'EDIT_DATASET_DIALOG';
export const IMPORT_DATASET_DIALOG = 'IMPORT_DATASET_DIALOG';
export const SAVE_DATASET_FILTER_DIALOG = 'SAVE_DATASET_FILTER_DIALOG';
export const DELETE_DATASET_DIALOG = 'DELETE_DATASET_DIALOG';

export const SET_DATASET_CHOICES = 'SET_DATASET_CHOICES';
export const RESTORE_VIEW = 'RESTORE_VIEW';

export const SET_DOT_PLOT_DATA = 'SET_DOT_PLOT_DATA';
export const SET_DOT_PLOT_SORT_ORDER = 'SET_DOT_PLOT_SORT_ORDER';
export const SET_EMBEDDING_DATA = 'SET_EMBEDDING_DATA';

export const SET_LOADING = 'SET_LOADING';

export const SET_LOADING_APP = 'LOADING_APP';

export const SET_NUMBER_OF_BINS_UI = 'SET_NUMBER_OF_BINS_UI';
export const SET_MARKER_SIZE_UI = 'SET_MARKER_SIZE_UI';
export const SET_MARKER_OPACITY_UI = 'SET_MARKER_OPACITY_UI';
export const SET_UNSELECTED_MARKER_OPACITY_UI = 'SET_UNSELECTED_MARKER_OPACITY_UI';


function getUser() {
    return function (dispatch, state) {
        fetch(API + '/user', {
            headers: {
                'Authorization': 'Bearer ' + getIdToken(),
                'Content-Type': 'application/json'
            }
        }).then(result => result.json()).then(user => dispatch(setUser(user)));
    };
}


export function initGapi() {
    return function (dispatch, state) {
        fetch(API + '/server').then(result => result.json()).then(serverInfo => {
            dispatch(setServerInfo(serverInfo));
            if (serverInfo.clientId == '') { // serving files locally, no login required
                dispatch(_setLoadingApp(false));
                dispatch(_setEmail(''));
                dispatch(listDatasets()).then(() => {
                    dispatch(_loadSavedView());
                });
            } else {
                let script = document.createElement('script');
                script.type = 'text/javascript';
                script.src = 'https://apis.google.com/js/api.js';
                script.onload = (e) => {
                    window.gapi.load('client:auth2', () => {
                        window.gapi.client.init({
                            clientId: serverInfo.clientId,
                            discoveryDocs: [],
                            scope: authScopes.join(' '),
                        }).then(() => {
                            dispatch(_setLoadingApp(false));
                            dispatch(initLogin(true));
                        });
                    });
                };
                document.getElementsByTagName('head')[0].appendChild(script);
            }
        }).catch(err => {
            handleError(dispatch, err);
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

        fetch(API + '/filter?id=' + filterId,
            {
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(response => response.json())
            .then(result => {
                result.id = filterId;
                result.value = JSON.parse(result.value);
                if (result.value.selected_points) {
                    dispatch(setDatasetFilter({}));
                    dispatch(handleFilterUpdated(result.value.selected_points));
                } else {
                    dispatch(setDatasetFilter(result.value));
                    dispatch(handleFilterUpdated());
                }
                dispatch(setSavedDatasetFilter(cloneDeep(result)));
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
        fetch(API + '/filter',
            {
                body: JSON.stringify({id: filterId}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            })
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
                let savedDatasetFilter = getState().savedDatasetFilter;
                if (savedDatasetFilter.id === filterId) {
                    dispatch(setSavedDatasetFilter(null));
                }
            }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to delete filter. Please try again.');
        });
    };
}

export function saveDatasetFilter(payload) {
    return function (dispatch, getState) {
        const state = getState();
        let filterValue;
        const datasetId = getState().dataset.id;

        if (Object.keys(state.datasetFilter).length > 0) {
            filterValue = state.datasetFilter;
        } else {
            filterValue = {selected_points: state.selection.points};
        }
        let savedDatasetFilter = state.savedDatasetFilter;
        let isEditFilter = savedDatasetFilter != null;
        let requestBody = {ds_id: datasetId};
        if (isEditFilter) {
            requestBody.id = savedDatasetFilter.id;
        }
        if (savedDatasetFilter == null) {
            savedDatasetFilter = {};
        }
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
        fetch(API + '/filter',
            {
                body: JSON.stringify(requestBody),
                method: isEditFilter ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(response => response.json())
            .then(result => {
                requestBody.id = result.id;
                let datasetFilters = getState().datasetFilters;
                if (isEditFilter) {
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
                dispatch(setSavedDatasetFilter(Object.assign({}, savedDatasetFilter)));
            }).finally(() => {
            dispatch(_setLoading(false));
            dispatch(setDialog(null));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to save filter. Please try again.');
        });
    };
}

export function removeDatasetFilter(filter) {
    return function (dispatch, getState) {
        if (filter == null) {
            // clear all
            dispatch(setDatasetFilter({}));
            dispatch(setSelection({}));
            dispatch(setSavedDatasetFilter(null));
        } else if (filter[0] === 'selection') {
            dispatch(setSelection({}));
        } else {
            let datasetFilter = getState().datasetFilter;
            delete datasetFilter[filter[0]];
            dispatch(setDatasetFilter(Object.assign({}, datasetFilter)));
        }
        dispatch(handleFilterUpdated());
    };
}

export function getDatasetFilterArray(datasetFilter) {
    let filters = [];
    for (let key in datasetFilter) {
        const value = datasetFilter[key];
        let f = null;
        if (window.Array.isArray(value)) {
            f = [key, 'in', value];
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

function getBinInfo(state) {
    if (state.dataset != null && state.dataset.summary != null && state.dataset.summary.embeddings != null) {
        // TODO handle more than one precomputed embedding
        const precomputedEmbedding = state.dataset.summary.embeddings[0];
        //json.embedding = precomputedEmbedding.basis;
        return {nbins: precomputedEmbedding.nbins, agg: precomputedEmbedding.agg};
    } else if (state.binValues) {
        return {nbins: state.numberOfBins, agg: state.binSummary};
    }
}

function getSliceJson(state, includeFilter = true) {
    const datasetId = state.dataset.id;
    const selectedEmbeddings = state.embeddings;
    let json = {
        id: datasetId,
        embedding: selectedEmbeddings.length > 0 ? selectedEmbeddings[0] : null,
    };
    let binInfo = getBinInfo(state);
    if (binInfo != null) {
        json.nbins = binInfo.nbins;
        json.agg = binInfo.agg;
    }

    if (includeFilter) {
        let filters = getDatasetFilterArray(state.datasetFilter);
        if (filters.length > 0) {
            const obs = state.dataset.obs;
            const obsCat = state.dataset.obsCat;
            for (let i = 0; i < filters.length; i++) {
                if (obsCat.indexOf(filters[i][0]) !== -1 || obs.indexOf(filters[i][0]) !== -1) {
                    filters[i][0] = 'obs/' + filters[i][0];
                }
            }
            json.filter = {filters: filters};
        } else {
            if (state.selection.points != null) {
                json.filter = {selected_points: state.selection.points};
            }
        }
    }
    return json;
}

export function downloadSelectedIds() {
    return function (dispatch, getState) {
        dispatch(_setLoading(true));
        let json = getSliceJson(getState());
        json['types'] = ['ids'];
        fetch(API + '/slice',
            {
                body: JSON.stringify(json),
                method: 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json()).then(result => {
            const blob = new Blob([result.join('\n')], {type: "text/plain;charset=utf-8"});
            saveAs(blob, "selection.txt");
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

function setSavedDatasetFilter(payload) {
    return {type: SET_SAVED_DATASET_FILTER, payload: payload};
}


/**
 *
 * @param selectedPoints Selected chart points if selection was updated or null.
 * @returns {Function}
 */
function handleFilterUpdated(selectedPoints) {
    return function (dispatch, getState) {

        // whenever filter is updated, we need to get selection statistics
        // if filter is not brush filter, we also need to get coordinates of selected cells
        const state = getState();
        let json = getSliceJson(state);
        if (json.filter != null && json.filter.selected_points != null) {
            json.filter = null; // do not use current brush selection
        }
        if (selectedPoints != null && selectedPoints.length > 0) {
            json.filter = {selected_points: selectedPoints};
        }
        const binInfo = getBinInfo(state);
        const traceBins = binInfo && state.embeddingData.length > 0 ? state.embeddingData[0].data[0].bins : null;
        json.types = ['selectionSummary', 'selectionCoordinates'];
        json.summary_measures = state.features;
        json.summary_dimensions = state.groupBy;

        if (json.filter == null) {
            dispatch(setSelection({}));
            return dispatch(setFeatureSummary({}));
        }
        dispatch(_setLoading(true));
        fetch(API + '/slice',
            {
                body: JSON.stringify(json),
                method: 'POST',
                headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json()).then(result => {
            const selectionSummary = result.selectionSummary;
            const selectionCoordinates = result.selectionCoordinates;
            if (selectionCoordinates != null) {
                const selectedIndicesOrBins = selectionCoordinates.indices_or_bins;
                let selectedPoints = selectedIndicesOrBins;
                if (binInfo != null) {
                    if (traceBins != null) {
                        selectedPoints = PlotUtil.convertBinsToPoints(traceBins, selectedIndicesOrBins);
                    } else {
                        console.log('Unable to convert bins to points');
                    }
                }

                dispatch(setSelection({
                    userPoints: selectedPoints,
                    points: selectedIndicesOrBins,
                    count: selectionCoordinates.count
                }));
            } else {
                dispatch(setSelection({}));
            }
            // userPoints are in chart space, points are in server space, count is total number of cells selected
            dispatch(setFeatureSummary(selectionSummary));
        }).catch(err => {
            handleError(dispatch, err);
        }).finally(() => dispatch(_setLoading(false)));
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
        dispatch(setSelection({}));
        dispatch(setDatasetFilter(Object.assign({}, datasetFilter)));
        if (update) {
            dispatch(handleFilterUpdated());
        }
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
            if (embeddingData[i].data[0].name === name) {
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
        dispatch(setSelection({}));
        dispatch(setDatasetFilter(datasetFilter));
        dispatch(handleFilterUpdated());
    };
}

export function handleSelectedPoints(payload) {
    return function (dispatch, getState) {
        let isEmpty = payload == null || payload.points.length === 0;
        let selectedpoints = isEmpty ? [] : payload.points[0].data.selectedpoints;
        if (!isEmpty && payload.points[0].data.bins != null) {
            selectedpoints = PlotUtil.convertPointsToBins(selectedpoints, payload.points[0].data.bins);
        }
        dispatch(setDatasetFilter({}));
        dispatch(handleFilterUpdated(selectedpoints));

    };
}

export function restoreView(payload) {
    return {type: RESTORE_VIEW, payload: payload};
}

export function setNumberOfBinsUI(payload) {
    return {type: SET_NUMBER_OF_BINS_UI, payload: payload};
}

export function setMarkerSize(payload) {
    return {type: SET_MARKER_SIZE, payload: payload};
}

export function setMarkerSizeUI(payload) {
    return {type: SET_MARKER_SIZE_UI, payload: payload};
}


export function setEmbeddingChartSize(payload) {
    return {type: SET_EMBEDDING_CHART_SIZE, payload: payload};
}


export function setUnselectedMarkerSize(payload) {
    return {type: SET_UNSELECTED_MARKER_SIZE, payload: payload};
}

export function setUnselectedMarkerSizeUI(payload) {
    return {type: SET_UNSELECTED_MARKER_SIZE_UI, payload: payload};
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
        const names = dataset.embeddings.map(e => e.name);
        let priority = {'X_fle': 1, 'X_umap': 2, 'X_tsne': 3, 'X_fitsne': 4};
        names.sort((a, b) => {
            a = priority[a] || Number.MAX_VALUE;
            b = priority[b] || Number.MAX_VALUE;
            return a - b;
        });

        if (names.length > 0) {
            dispatch(setSelectedEmbedding([names[0]]));
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
            if (savedView.colorScheme != null) {
                let interp = getInterpolator(savedView.colorScheme);
                if (interp != null) {
                    savedView.colorScheme = {
                        name: savedView.colorScheme,
                        value: interp
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
            dispatch(setDataset(savedView.dataset, false, false))
                .then(() => dispatch(setDatasetFilter(savedView.datasetFilter)))
                .then(() => dispatch(restoreView(savedView)))
                .then(() => dispatch(_updateCharts({dotplot: true, clear: true})))
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

// export function getAccessToken() {
//     return window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().access_token;
// }

export function getIdToken() {
    return typeof window.gapi !== 'undefined' ? window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().id_token : '';
}

export function saveDataset(payload) {
    return function (dispatch, getState) {
        const readers = payload.readers;
        const name = payload.name;
        const url = payload.url;

        if (name == '' || url === '') {
            return;
        }
        // let bucket = url.substring('gs://'.length);
        // let slash = bucket.indexOf('/');
        // let object = encodeURIComponent(bucket.substring(slash + 1));
        // bucket = encodeURIComponent(bucket.substring(0, slash));
        let isEdit = payload.dataset != null;
        dispatch(_setLoading(true));

        const serverEmail = getState().serverInfo.email;
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
        }).then(() => fetch(API + '/dataset',
            {
                body: JSON.stringify(
                    {
                        id: payload.dataset != null ? payload.dataset.id : null,
                        url: url,
                        name: name,
                        readers: readers
                    }),
                method: isEdit ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            })).then(importDatasetResponse => importDatasetResponse.json()).then(importDatasetResult => {
            if (isEdit) {
                dispatch(updateDataset({name: name, id: importDatasetResult.id, owner: true}));
                dispatch(setMessage('Dataset updated'));
            } else {
                dispatch(_addDataset({name: name, id: importDatasetResult.id, owner: true}));
                dispatch(setMessage(updatePermissions ? 'Please ensure that ' + serverEmail + ' is a "Storage Object Viewer"' : 'Dataset created'));
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
        fetch(API + '/dataset',
            {
                body: JSON.stringify(
                    {id: payload.dataset.id}),
                method: 'DELETE',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(() => {
            dispatch(_deleteDataset({id: payload.dataset.id}));
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

function _setEmbeddingData(payload) {
    return {type: SET_EMBEDDING_DATA, payload: payload};
}

export function setDialog(payload) {
    return {type: SET_DIALOG, payload: payload};
}

function _setDatasetChoices(payload) {
    return {type: SET_DATASET_CHOICES, payload: payload};
}

function _setLoading(payload) {
    return {type: SET_LOADING, payload: payload};
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

function _setEmail(payload) {
    return {type: SET_EMAIL, payload: payload};
}


export function setFeatures(payload) {
    const value = payload;
    return function (dispatch, getState) {
        const state = getState();
        const obsCat = state.dataset.obsCat;
        const priorFeatures = state.features;
        const priorGroupBy = state.groupBy;
        const features = [];
        const groupBy = [];
        value.forEach(val => {
            if (obsCat.indexOf(val) !== -1) {
                groupBy.push(val);
            } else {
                features.push(val);
            }
        });
        if (features.length !== priorFeatures.length) {
            dispatch(_setFeatures(features));
        }
        if (groupBy.length !== priorGroupBy.length) {
            dispatch(_setGroupBy(groupBy));
        }
    };
}


function _setFeatures(payload) {

    return function (dispatch, getState) {
        let prior = getState().features;
        dispatch({type: SET_FEATURES, payload: payload});
        dispatch(_updateCharts({dotplot: true}, err => {
            dispatch({type: SET_FEATURES, payload: prior});
        }));
    };
}

function _setGroupBy(payload) {
    return function (dispatch, getState) {
        let prior = getState().groupBy; // in case of error, restore
        // let datasetFilter = getState().datasetFilter;
        // let categoricalFilterChanged = false;
        // for (let key in datasetFilter) {
        //     if (payload.indexOf(key) === -1) {
        //         delete datasetFilter[key];
        //         categoricalFilterChanged = true;
        //     }
        // }
        // if (categoricalFilterChanged) {
        //     dispatch({type: SET_DATASET_FILTER, payload: datasetFilter});
        // }
        dispatch({type: SET_GROUP_BY, payload: payload}); // updated choices
        dispatch(_updateCharts({dotplot: true}, err => {
            dispatch({type: SET_GROUP_BY, payload: prior});
        }));
    };
}

export function setSelectedEmbedding(payload) {
    return function (dispatch, getState) {
        let prior = getState().embeddings;
        dispatch({type: SET_SELECTED_EMBEDDING, payload: payload});
        dispatch(_updateCharts({clear: true}, err => {
            dispatch({type: SET_SELECTED_EMBEDDING, payload: prior});
        }));
    };
}


export function setNumberOfBins(payload) {
    return function (dispatch, getState) {
        if (getState().numberOfBins !== payload) {
            let prior = getState().numberOfBins;
            dispatch({type: SET_NUMBER_OF_BINS, payload: payload});
            dispatch(_updateCharts({clear: true}, err => {
                dispatch({type: SET_NUMBER_OF_BINS, payload: prior});
            }));
        }
    };
}

export function setBinSummary(payload) {
    return function (dispatch, getState) {
        let prior = getState().binSummary;
        dispatch({type: SET_BIN_SUMMARY, payload: payload});

        dispatch(_updateCharts({clear: true}, err => {
            dispatch({type: SET_BIN_SUMMARY, payload: prior});
        }));
    };
}

export function setBinValues(payload) {
    return function (dispatch, getState) {
        let prior = getState().binValues;
        dispatch({type: SET_BIN_VALUES, payload: payload});
        dispatch(_updateCharts({clear: true}, err => {
            dispatch({type: SET_BIN_VALUES, payload: prior});
        }));

    };
}


function setDatasetFilters(payload) {
    return {type: SET_DATASET_FILTERS, payload: payload};
}


export function setDataset(id, loadDefaultView = true, setLoading = true) {
    return function (dispatch, getState) {
        const datasetChoices = getState().datasetChoices;
        let choice = null;
        for (let i = 0; i < datasetChoices.length; i++) {
            if (datasetChoices[i].id === id) {
                choice = datasetChoices[i];
                break;
            }
        }
        if (choice == null) {
            dispatch(setMessage('Unable to find dataset'));
            return Promise.reject('Unable to find dataset');
        }
        // force re-render selected dataset dropdown
        dispatch(_setDataset({
            owner: choice.owner,
            name: choice.name,
            id: id,
            embeddings: [],
            features: [],
            obs: [],
            obsCat: []
        }));
        if (setLoading) {
            dispatch(_setLoading(true));
        }

        const filtersPromise = fetch(API + '/filters?id=' + id,
            {
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json())
            .then(results => {
                dispatch(setDatasetFilters(results));
            });
        const schemaPromise = fetch(API + '/schema?id=' + id, {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        }).then(result => {
            let newDataset = result;
            newDataset.owner = choice.owner;
            newDataset.name = choice.name;
            newDataset.features = result.var;
            newDataset.id = id;
            dispatch(_setDataset(newDataset));
            if (loadDefaultView) {
                dispatch(loadDefaultDatasetEmbedding());
            }
        });
        return Promise.all([schemaPromise, filtersPromise]).finally(() => {
            if (setLoading) {
                dispatch(_setLoading(false));
            }
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve dataset. Please try again.');
        });
    };
}


// dot plot depends on features, groupBy
// embedding chart depends on groupBy, features, embedding (also markerSize, colorScale which don't require API
// call). Also need to updated embedding when selection changes.
/**
 *
 * @param sliceOptions.dotplot
 * @param sliceOptions.clear Clear cached data, such as when # of bins changes
 * @returns {function(*, *): Promise<any | never>}
 */
function _updateCharts(sliceOptions, onError) {

    return function (dispatch, getState) {
        const state = getState();
        if (state.dataset == null) {
            return;
        }
        const continuousFeatures = state.features;
        const categoricalFeatures = state.groupBy;
        dispatch(_setLoading(true));
        const obs = state.dataset.obs;
        const selectedEmbeddings = state.embeddings;
        if (sliceOptions.dotplot && continuousFeatures.length === 0) {
            sliceOptions.dotplot = false;
        }

        let interpolator = state.interpolator;
        let rgbScale = scaleLinear().domain([0, 255]).range([0, 1]);
        const markerSize = state.markerSize;
        const markerOpacity = state.markerOpacity;
        const unselectedMarkerSize = state.unselectedMarkerSize;
        const unselectedMarkerOpacity = state.unselectedMarkerOpacity;
        const obsCat = state.dataset.obsCat;
        const dotPlotData = state.dotPlotData;
        let embeddingData = state.embeddingData;
        const embeddingChartSize = state.embeddingChartSize;
        if (sliceOptions.clear) {
            embeddingData = [];
        }
        let cachedFeatureNames = {};
        let currentFeatures = continuousFeatures.concat(categoricalFeatures);
        if (currentFeatures.length === 0) {
            currentFeatures = ['__count'];
        }
        // reuse existing traces
        embeddingData.forEach(traceInfo => {
            cachedFeatureNames[traceInfo.name] = true;
            let active = currentFeatures.indexOf(traceInfo.name) !== -1;
            if (active) {
                traceInfo.date = new Date();
                traceInfo.layout = PlotUtil.createEmbeddingLayout(
                    {
                        size: embeddingChartSize,
                        is3d: selectedEmbeddings[0].endsWith('3d'),
                        title: traceInfo.name
                    });
            }
            traceInfo.active = active;
        });
        let embeddingJson = getSliceJson(state, false);
        embeddingJson.types = [];
        embeddingJson.embedding_measures = [];
        embeddingJson.embedding_dimensions = [];
        let dotPlotJson = Object.assign({}, embeddingJson);
        dotPlotJson.types = [];
        let dotPlotKeys = {};
        if (sliceOptions.dotplot) {
            dotPlotJson.dotplot_measures = [];
            dotPlotJson.dotplot_dimensions = [];
            let cachedDotPlotKeys = {};
            dotPlotData.forEach(dotPlotDataItem => {
                dotPlotDataItem.values.forEach(datum => {
                    cachedDotPlotKeys[dotPlotDataItem.name + '-' + datum.name] = true;
                });
            });

            categoricalFeatures.forEach(category => {
                let added = false;
                continuousFeatures.forEach(feature => {
                    let isObs = obs.indexOf(feature) !== -1;
                    if (!isObs) {
                        let key = category + '-' + feature;
                        dotPlotKeys[key] = true;
                        if (cachedDotPlotKeys[key] == null) {
                            dotPlotJson.dotplot_measures.push(feature);
                            added = true;
                        }
                    }
                });
                if (added) {
                    dotPlotJson.dotplot_dimensions.push(category);
                }
            });
        }
        categoricalFeatures.forEach(feature => {
            let isCached = cachedFeatureNames[feature] != null;
            if (!isCached) {
                embeddingJson.embedding_dimensions.push(feature);
            }
        });

        continuousFeatures.forEach(feature => {
            let isObs = obs.indexOf(feature) !== -1;
            let isCached = cachedFeatureNames[feature] != null;
            if (!isCached) {
                embeddingJson.embedding_measures.push(isObs ? 'obs.' + feature : feature);
            }
        });


        if (sliceOptions.dotplot && dotPlotJson.dotplot_dimensions.length > 0 && dotPlotJson.dotplot_measures.length > 0) {
            dotPlotJson.types = ['dotplot'];
        }
        const numEmbeddingFeaturesRequested = embeddingJson.embedding_measures.length + embeddingJson.embedding_dimensions.length;
        if (numEmbeddingFeaturesRequested > 0 || cachedFeatureNames['__count'] == null) {
            embeddingJson.types = ['embedding', 'summary'];
            if (numEmbeddingFeaturesRequested === 0) {
                embeddingJson.embedding_measures = ['__count'];
            }
            embeddingJson.summary_measures = embeddingJson.embedding_measures;
            embeddingJson.summary_dimensions = embeddingJson.embedding_dimensions;
        }
        let embeddingPromise;
        let dotPlotPromise;
        let selectionSummaryPromise;
        if (embeddingJson.types.length === 0) {
            embeddingPromise = Promise.resolve({});
            selectionSummaryPromise = Promise.resolve({});
        } else {
            let selectionSummaryJson = getSliceJson(state);
            selectionSummaryJson.types = ['selectionSummary'];
            selectionSummaryJson.summary_measures = embeddingJson.embedding_measures;
            selectionSummaryJson.summary_dimensions = embeddingJson.embedding_dimensions;
            if ((selectionSummaryJson.summary_measures.length == 1 && selectionSummaryJson.summary_measures[0] === '__count') || selectionSummaryJson.filter == null) {
                selectionSummaryPromise = Promise.resolve({});
            } else {

                selectionSummaryPromise = fetch(API + '/slice',
                    {
                        body: JSON.stringify(selectionSummaryJson),
                        method: 'POST',
                        headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
                    }).then(r => r.json());
            }
            embeddingPromise = fetch(API + '/slice',
                {
                    body: JSON.stringify(embeddingJson),
                    method: 'POST',
                    headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
                }).then(r => r.json());
        }
        if (dotPlotJson.types.length === 0) {
            dotPlotPromise = Promise.resolve({});
        } else {
            dotPlotPromise = fetch(API + '/slice',
                {
                    body: JSON.stringify(dotPlotJson),
                    method: 'POST',
                    headers: {'Content-Type': 'application/json', 'Authorization': 'Bearer ' + getIdToken()},
                }).then(r => r.json());
        }

        selectionSummaryPromise.then(result => {
            if (result.selectionSummary) {
                let selectionSummary = Object.assign(getState().featureSummary, result.selectionSummary);
                dispatch(setFeatureSummary(selectionSummary));
            }
        });

        dotPlotPromise.then(dotPlotResult => {
            if (dotPlotResult.dotplot != null) {
                let newDotplotData = dotPlotResult.dotplot;
                // merge with existing data and set active flags
                newDotplotData.forEach(newDotplotDataItem => {
                    let existingIndex = -1;
                    for (let i = 0; i < dotPlotData.length; i++) {
                        if (dotPlotData[i].name === newDotplotDataItem.name) {
                            existingIndex = i;
                            break;
                        }
                    }
                    if (existingIndex === -1) {
                        dotPlotData.push(newDotplotDataItem);
                    } else {
                        newDotplotDataItem.values.forEach(value => {
                            dotPlotData[existingIndex].values.push(value);
                        });
                    }
                });
            }
            dotPlotData.forEach((dotPlotDataItem, dotPlotDataItemIndex) => {
                let active = categoricalFeatures.indexOf(dotPlotDataItem.name) !== -1;
                dotPlotDataItem.active = active;
                if (active) {
                    let oneActiveFeature = false;
                    dotPlotDataItem.values.forEach(datum => {
                        datum.active = dotPlotKeys[dotPlotDataItem.name + '-' + datum.name];
                        if (datum.active) {
                            oneActiveFeature = true;
                        }
                    });
                    dotPlotDataItem.active = oneActiveFeature;
                    if (dotPlotDataItem.active) {
                        dotPlotData[dotPlotDataItemIndex] = Object.assign({}, dotPlotDataItem);
                        dotPlotDataItem.values.sort((a, b) => {
                            a = a.name;
                            b = b.name;
                            return (a < b ? -1 : (a === b ? 0 : 1));
                        });
                    }
                }
            });
            dotPlotData.sort((a, b) => {
                a = a.name;
                b = b.name;
                return (a < b ? -1 : (a === b ? 0 : 1));
            });
            dispatch(_setDotPlotData(dotPlotData.slice()));
        });
        embeddingPromise.then(sliceResult => {

                const embeddingResult = sliceResult.embedding || {};
                const embeddingBins = embeddingResult.bins;
                const embeddingValues = embeddingResult.values;
                const newSummary = sliceResult.summary || {};
                const globalFeatureSummary = state.globalFeatureSummary;
                for (let key in newSummary) {
                    globalFeatureSummary[key] = newSummary[key];
                }
                const coordinates = embeddingResult.coordinates;
                const is3d = selectedEmbeddings[0].endsWith('3d');
                const userPoints = state.selection.userPoints;

                // add new embedding values
                for (let name in embeddingValues) {
                    let traceSummary = globalFeatureSummary[name];
                    let x = coordinates[selectedEmbeddings[0] + '_1'];
                    let y = coordinates[selectedEmbeddings[0] + '_2'];
                    let z = coordinates[selectedEmbeddings[0] + '_3'];
                    let values = embeddingValues[name];
                    let purity = null;
                    if (isPlainObject(values)) {
                        values = values.mode;
                        purity = values.purity;
                    }

                    let isCategorical = name !== '__count' && obsCat.indexOf(name) !== -1;
                    let colorScale = null;

                    if (!isCategorical) {
                        if (traceSummary == null) { // __count is not currently precomputed
                            console.log('Computing range for ' + name);
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
                        colorScale = scaleSequential(interpolator.value).domain([traceSummary.min, traceSummary.max]);
                        colorScale.summary = traceSummary;
                    } else {
                        let traceUniqueValues = traceSummary.categories;
                        colorScale = scaleOrdinal(
                            traceUniqueValues.length <= 10 ? schemeCategory10 : (traceUniqueValues.length <= 20 ? CATEGORY_20B : CATEGORY_20B.concat(
                                CATEGORY_20C))).domain(traceUniqueValues);
                        colorScale.summary = traceSummary;
                    }

                    let colors = [];
                    for (let i = 0, n = values.length; i < n; i++) {
                        let rgb = color(colorScale(values[i]));
                        colors.push([rgbScale(rgb.r), rgbScale(rgb.g), rgbScale(rgb.b)]);
                    }

                    let chartData = {
                        hoverinfo: 'text',
                        showlegend: false,
                        name: name,
                        mode: 'markers',
                        type: is3d ? 'scatter3d' : 'scattergl',
                        hoveron: 'points',
                        x: x,
                        y: y,
                        bins: embeddingBins,
                        marker: {
                            size: markerSize,
                            color: colors,
                            opacity: markerOpacity,
                            showscale: false,
                        },
                        unselected: {marker: {opacity: unselectedMarkerOpacity, size: unselectedMarkerSize}},
                        values: values,
                        purity: purity,
                        text: values,
                    };
                    if (is3d) {
                        chartData.z = z;
                    }

                    chartData = [chartData];
                    embeddingData.push(
                        {
                            date: new Date(),
                            active: true,
                            colorScale: colorScale,
                            continuous: !isCategorical,
                            data: chartData,
                            name: name,
                            layout: PlotUtil.createEmbeddingLayout(
                                {size: embeddingChartSize, is3d: is3d, title: name}),
                            isCategorical: isCategorical,
                        });
                }
                embeddingData.sort((a, b) => {
                    a = a.name.toLowerCase();
                    b = b.name.toLowerCase();
                    return a < b ? -1 : 1;
                });
                embeddingData.forEach(item => {
                    item.data.forEach(trace => {
                        trace.selectedpoints = userPoints;
                    });
                });
                dispatch(setGlobalFeatureSummary(globalFeatureSummary));
                dispatch(_setEmbeddingData(embeddingData.slice(0)));
            },
        );
        return Promise.all([embeddingPromise, dotPlotPromise, selectionSummaryPromise]).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve data. Please try again.');
            if (onError) {
                onError(err);
            }
        });

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
        return fetch(API + '/datasets', {headers: {'Authorization': 'Bearer ' + getIdToken()}})
            .then(response => {
                return response.json();
            })
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



import {color} from 'd3-color';
import {scaleLinear, scaleOrdinal, scaleSequential} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import CustomError from '../CustomError';
import PlotUtil, {getInterpolator} from '../PlotUtil';


export const API = '/api';
const authScopes = [
    'email',
    'profile',
    'https://www.googleapis.com/auth/userinfo.profile',
    'https://www.googleapis.com/auth/contacts.readonly',
    'https://www.googleapis.com/auth/devstorage.full_control',
];

export const SET_UNSELECTED_MARKER_SIZE = 'SET_UNSELECTED_MARKER_SIZE';
export const SET_UNSELECTED_MARKER_SIZE_UI = 'SET_UNSELECTED_MARKER_SIZE_UI';
export const SET_SERVER_INFO = "SET_SERVER_INFO";
export const SET_LEGEND_VISIBILITY = 'SET_LEGEND_VISIBILITY';
export const ADD_DATASET = 'ADD_DATASET';
export const DELETE_DATASET = 'DELETE_DATASET';
export const UPDATE_DATASET = 'UPDATE_DATASET';

export const SET_MARKER_SIZE = 'SET_MARKER_SIZE';
export const SET_MARKER_OPACITY = 'SET_MARKER_OPACITY';

export const SET_EMBEDDING_CHART_SIZE = "SET_EMBEDDING_CHART_SIZE";

export const SET_UNSELECTED_MARKER_OPACITY = 'SET_UNSELECTED_MARKER_OPACITY';

// update chart

export const SET_SELECTED_VALUE_COUNTS = 'SET_SELECTED_VALUE_COUNTS';
export const SET_FEATURES = 'SET_FEATURES';
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
export const DELETE_DATASET_DIALOG = 'DELETE_DATASET_DIALOG';

export const SET_DATASET_CHOICES = 'SET_DATASET_CHOICES';
export const RESTORE_VIEW = 'RESTORE_VIEW';

export const SET_DOT_PLOT_DATA = 'SET_DOT_PLOT_DATA';
export const SET_EMBEDDING_DATA = 'SET_EMBEDDING_DATA';

export const SET_LOADING = 'SET_LOADING';

export const SET_LOADING_APP = 'LOADING_APP';

export const SET_NUMBER_OF_BINS_UI = 'SET_NUMBER_OF_BINS_UI';
export const SET_MARKER_SIZE_UI = 'SET_MARKER_SIZE_UI';
export const SET_MARKER_OPACITY_UI = 'SET_MARKER_OPACITY_UI';
export const SET_UNSELECTED_MARKER_OPACITY_UI = 'SET_UNSELECTED_MARKER_OPACITY_UI';

const TWENTY_COLORS = [
    '#1f77b4', '#aec7e8', '#ff7f0e',
    '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd',
    '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
    '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5'];

const CATEGORY_20B = [
    '#393b79', '#5254a3', '#6b6ecf',
    '#9c9ede', '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31',
    '#bd9e39', '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b',
    '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6'];
const CATEGORY_20C = [
    '#3182bd', '#6baed6', '#9ecae1',
    '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2', '#31a354',
    '#74c476', '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc',
    '#dadaeb', '#636363', '#969696', '#bdbdbd', '#d9d9d9'];

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
            if (serverInfo.clientId == '') { // serving files locally, no login required
                dispatch(_setLoadingApp(false));
                dispatch(_setEmail(''));
                dispatch(listDatasets()).then(() => {
                    dispatch(_loadSavedView());
                    dispatch(loadDefaultDataset());
                });
            } else {
                let script = document.createElement('script');
                script.type = 'text/javascript';
                script.src = 'https://apis.google.com/js/api.js';
                dispatch(setServerInfo(serverInfo));
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


function setSelectedValueCounts(payload) {
    return {type: SET_SELECTED_VALUE_COUNTS, payload: payload};
}

function selectedValueCountsPromise(selectedpoints, legendVisibility, dispatch, getState) {
    const datasetId = getState().dataset.id;
    const continuousFeatures = getState().features;
    const categoricalFeatures = getState().groupBy;
    const selectedEmbeddings = getState().embeddings;

    let json = {
        id: datasetId,
        embedding: selectedEmbeddings[0],
        keys: continuousFeatures.concat(categoricalFeatures)
    };
    if (getState().binValues) {
        json.nbins = getState().numberOfBins;
    }
    if (selectedpoints.length > 0) {
        json.p = selectedpoints;
    }
    if (Object.keys(legendVisibility).length > 0) {
        json.c = legendVisibility;
    }
    let p;
    if (selectedEmbeddings.length === 0 || (selectedpoints.length === 0 && json.c == null)) {
        p = Promise.resolve({});
    } else {
        p = fetch(API + '/selected_value_counts',
            {
                body: JSON.stringify(json),
                method: 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            }).then(result => result.json());
    }
    p.then(selectedValueCounts => dispatch(setSelectedValueCounts(selectedValueCounts)));
    return p;
}


function setLegendVisibility(payload) {
    return {type: SET_LEGEND_VISIBILITY, payload: payload};
}

export function handleLegendClick(payload) {
    return function (dispatch, getState) {
        let name = payload.name;
        let value = payload.value;
        // legend visibility, true means hidden : {'leiden':{'1':true, '5':true}}
        // keep track of hidden values
        let legendVisibility = getState().legendVisibility;
        let values = legendVisibility[name];
        if (values == null) {
            values = [];
            legendVisibility[name] = values;
        }
        let index = values.indexOf(value);
        if (index !== -1) {
            values.splice(index, 1);
            if (values.length === 0) {
                delete legendVisibility[name];
            }
        } else {
            values.push(value);
        }

        dispatch(setLegendVisibility(legendVisibility));
        selectedValueCountsPromise([], legendVisibility, dispatch, getState);
    };
}

export function handleSelectedPoints(payload) {
    // we intentionally don't persist selected points
    let selectedpoints = payload == null ? [] : payload.points[0].data.selectedpoints;
    return function (dispatch, getState) {
        selectedValueCountsPromise(selectedpoints, getState().legendVisibility, dispatch, getState);
    };
    // TODO linked brushing across embeddings


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
        let q = window.location.search.substring(3);
        if (q.length === 0) {
            return;
        }
        let savedView;
        try {
            savedView = JSON.parse(window.decodeURIComponent(q));
        } catch (err) {
            return dispatch(setMessage('Unable to restore saved view.'));
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
            dispatch(setDataset(savedView.dataset))
                .then(() => dispatch(restoreView(savedView)))
                .then(() => dispatch(_updateEmbedding({
                    embedding: true,
                    dotPlot: true,
                    clear: true
                }))).catch(err => {
                console.log(err);
                dispatch(setMessage('Unable to restore saved view.'));
            });

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
        let bucket = url.substring('gs://'.length);
        let slash = bucket.indexOf('/');
        // let object = encodeURIComponent(bucket.substring(slash + 1));
        // bucket = encodeURIComponent(bucket.substring(0, slash));
        let isEdit = payload.dataset != null;
        dispatch(_setLoading(true));
        let updateDatasetPermissionPromise = null;
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
        updateDatasetPermissionPromise = Promise.resolve({ok: true});
        // }
        updateDatasetPermissionPromise.then(permissionsResponse => {

            // if (!permissionsResponse.ok) {
            //     dispatch(setMessage('Unable to set dataset read permissions. Please ensure that you are the dataset owner or manually add ' + serverEmail + ' as a reader.'));
            // }
        }).then(() => fetch(API + '/dataset',
            {
                body: JSON.stringify(
                    {id: payload.dataset != null ? payload.dataset.id : null, url: url, name: name, readers: readers}),
                method: isEdit ? 'PUT' : 'POST',
                headers: {'Authorization': 'Bearer ' + getIdToken()},
            })).then(importDatasetResponse => importDatasetResponse.json()).then(importDatasetResult => {
            if (isEdit) {
                dispatch(updateDataset({name: name, id: importDatasetResult.id, owner: true}));
                dispatch(setMessage(updatePermissions ? 'Please ensure that ' + serverEmail + ' is a "Storage Object Viewer"' : 'Dataset updated'));
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
        dispatch(_updateEmbedding({embedding: true, dotPlot: true}, err => {
            dispatch({type: SET_FEATURES, payload: prior});
        }));
    };
}

function _setGroupBy(payload) {
    return function (dispatch, getState) {
        let prior = getState().groupBy; // in case of error, restore
        let legendVisibility = getState().legendVisibility;
        let legendVisibilityChanged = false;
        for (let key in legendVisibility) {
            if (payload.indexOf(key) === -1) {
                delete legendVisibility[key];
                legendVisibilityChanged = true;
            }
        }
        if (legendVisibilityChanged) {
            dispatch({type: SET_LEGEND_VISIBILITY, payload: legendVisibility});
        }
        dispatch({type: SET_GROUP_BY, payload: payload}); // updated choices
        dispatch(_updateEmbedding({embedding: true, dotPlot: true}, err => {
            dispatch({type: SET_GROUP_BY, payload: prior});
        }));
    };
}

export function setSelectedEmbedding(payload) {
    return function (dispatch, getState) {
        let prior = getState().embeddings;
        dispatch({type: SET_SELECTED_EMBEDDING, payload: payload});
        dispatch(_updateEmbedding({embedding: true, clear: true}, err => {
            dispatch({type: SET_SELECTED_EMBEDDING, payload: prior});
        }));
    };
}


export function setNumberOfBins(payload) {
    return function (dispatch, getState) {
        if (getState().numberOfBins !== payload) {
            let prior = getState().numberOfBins;
            dispatch({type: SET_NUMBER_OF_BINS, payload: payload});
            dispatch(_updateEmbedding({embedding: true, clear: true}, err => {
                dispatch({type: SET_NUMBER_OF_BINS, payload: prior});
            }));
        }
    };
}

export function setBinSummary(payload) {
    return function (dispatch, getState) {
        let prior = getState().binSummary;
        dispatch({type: SET_BIN_SUMMARY, payload: payload});

        dispatch(_updateEmbedding({embedding: true, clear: true}, err => {
            dispatch({type: SET_BIN_SUMMARY, payload: prior});
        }));
    };
}

export function setBinValues(payload) {

    return function (dispatch, getState) {
        let prior = getState().binValues;
        dispatch({type: SET_BIN_VALUES, payload: payload});
        dispatch(_updateEmbedding({embedding: true, clear: true}, err => {
            dispatch({type: SET_BIN_VALUES, payload: prior});
        }));

    };
}

export function setDataset(id) {
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
            return;
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
        dispatch(_setLoading(true));
        let url = [API + '/schema?id=' + id];

        return fetch(url.join(''), {headers: {'Authorization': 'Bearer ' + getIdToken()}}).then(response => {
            return response.json();
        }).then(result => {
            let newDataset;
            if (result.version === '1') {
                newDataset = {
                    owner: choice.owner,
                    name: choice.name,
                    id: id,
                    embeddings: result.embeddings,
                    features: result.var,
                    obsCat: result.obs_cat,
                    nObs: result.n_obs,
                    obs: result.obs
                };
            }
            dispatch(_setDataset(newDataset));
            dispatch(loadDefaultDatasetEmbedding());
            return newDataset;


        }).finally(() => {
            dispatch(_setLoading(false));
        }).catch(err => {
            handleError(dispatch, err, 'Unable to retrieve dataset. Please try again.');
            // dispatch(_setDataset(null));
        });
    };
}

// dot plot depends on features, groupBy
// embedding chart depends on groupBy, features, embedding (also markerSize, colorScale which don't require API
// call)
/**
 *
 * @param options.embedding
 * @param options.dotPlot
 * @param options.clear Clear cached data
 * @returns {function(*, *): Promise<any | never>}
 */
function _updateEmbedding(options, onError) {

    return function (dispatch, getState) {

        const continuousFeatures = getState().features;
        const categoricalFeatures = getState().groupBy;

        if (getState().dataset == null) {
            return;
        }
        dispatch(_setLoading(true));
        const datasetId = getState().dataset.id;
        const obs = getState().dataset.obs;
        const selectedEmbeddings = getState().embeddings;
        if (options.embedding && selectedEmbeddings.length === 0) {
            options.embedding = false;
        }
        if (options.dotPlot && continuousFeatures.length === 0) {
            options.dotPlot = false;
        }

        let embeddingData = getState().embeddingData;
        const embeddingChartSize = getState().embeddingChartSize;
        if (options.clear) {
            embeddingData = [];
        }
        let cachedFeatureNames = {};
        const currentFeatures = continuousFeatures.concat(categoricalFeatures);

        let nRequestedFeatures = 0;
        const markerSize = getState().markerSize;
        const markerOpacity = getState().markerOpacity;
        const unselectedMarkerSize = getState().unselectedMarkerSize;

        const unselectedMarkerOpacity = getState().unselectedMarkerOpacity;

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

        let embeddingUrl = [API + '/slice?id=' + datasetId];
        let dotPlotUrl = [API + '/slice?id=' + datasetId];
        // TODO only get new features for dot plot
        // only get new features for embedding plot
        let ngenes = 0;
        continuousFeatures.forEach(feature => {
            if (obs.indexOf(feature) === -1) {
                dotPlotUrl.push('&key=' + feature);
                ngenes++;
            }
            let isCached = cachedFeatureNames[feature] != null;
            if (!isCached) {
                nRequestedFeatures++;
                embeddingUrl.push('&key=' + feature);
            }
        });

        if ((options.embedding || options.dotPlot)) {
            categoricalFeatures.forEach(feature => {
                dotPlotUrl.push('&key=' + feature);
                let isCached = cachedFeatureNames[feature] != null;
                if (!isCached) {
                    nRequestedFeatures++;
                    embeddingUrl.push('&key=' + feature);
                }
            });
        }

        let binValues = getState().numberOfBins > 0 && getState().binValues;
        if (options.embedding && binValues) {
            embeddingUrl.push('&nbins=' + getState().numberOfBins);
            embeddingUrl.push('&reduce_function=' + getState().binSummary);
        }

        if (options.embedding && selectedEmbeddings.length > 0) {
            selectedEmbeddings.forEach(selectedEmbedding => embeddingUrl.push('&embedding=' + selectedEmbedding));
        }

        let interpolator = getState().interpolator;
        let selectedValueCounts = getState().selectedValueCounts;
        let promises = [];
        if (options.dotPlot) {
            let dotPlotPromise = null;
            if (ngenes > 0 && categoricalFeatures.length > 0) {
                dotPlotPromise = fetch(dotPlotUrl.join(''),
                    {headers: {'Authorization': 'Bearer ' + getIdToken()}})
                    .then(response => {
                        return response.json();
                    });
            } else {
                dotPlotPromise = Promise.resolve({stats: {}});
            }
            promises.push(dotPlotPromise);
            dotPlotPromise.then(result => { // TODO cache
                dispatch(_setDotPlotData(result.dotplot || null));

            });
        }

        if (options.embedding) {
            let valueCountsPromise = selectedValueCountsPromise([], getState().legendVisibility, dispatch, getState);
            let embeddingPromise;
            if (nRequestedFeatures === 0 && embeddingData.filter(d => d.active).length > 0) {
                embeddingPromise = Promise.resolve({embedding: {values: {}}});
            } else {
                embeddingPromise = fetch(embeddingUrl.join(''), {headers: {'Authorization': 'Bearer ' + getIdToken()}})
                    .then(response => {
                        return response.json();
                    });
            }
            promises.push(embeddingPromise);
            promises.push(valueCountsPromise);
            let rgbScale = scaleLinear().domain([0, 255]).range([0, 1]);
            embeddingPromise.then(embeddingResult => {
                    const embeddingValues = embeddingResult.embedding.values;
                    const categories = embeddingResult.embedding.categories;
                    let coordinates = {};
                    const is3d = selectedEmbeddings[0].endsWith('3d');

                    // TODO > 1 embedding
                    coordinates[selectedEmbeddings[0] + '_1'] = new Float32Array(embeddingValues[selectedEmbeddings[0] + '_1']);
                    coordinates[selectedEmbeddings[0] + '_2'] = new Float32Array(embeddingValues[selectedEmbeddings[0] + '_2']);
                    if (is3d) {
                        coordinates[selectedEmbeddings[0] + '_3'] = new Float32Array(embeddingValues[selectedEmbeddings[0] + '_3']);
                    }

                    for (let name in embeddingValues) {
                        if (name in coordinates) { // skip coordinate columns and cell ids
                            continue;
                        }

                        let x = coordinates[selectedEmbeddings[0] + '_1'];
                        let y = coordinates[selectedEmbeddings[0] + '_2'];
                        let z = coordinates[selectedEmbeddings[0] + '_3'];
                        let values = embeddingValues[name];
                        let isCategorical = categories[name] != null;
                        let colorScale = null;
                        let min = null;
                        let max = null;
                        if (!isCategorical) {
                            min = Number.MAX_VALUE;
                            max = -Number.MAX_VALUE;
                            for (let i = 0, n = values.length; i < n; i++) {
                                let value = values[i];
                                min = value < min ? value : min;
                                max = value > max ? value : max;
                            }
                            colorScale = scaleSequential(interpolator.value).domain([min, max]);
                            colorScale.valueCounts = embeddingResult.embedding.obs[name];
                        } else {
                            let traceUniqueValues = categories[name].values;
                            colorScale = scaleOrdinal(
                                traceUniqueValues.length <= 10 ? schemeCategory10 : (traceUniqueValues.length <= 20 ? CATEGORY_20B : CATEGORY_20B.concat(
                                    CATEGORY_20C))).domain(traceUniqueValues);
                            colorScale.valueCounts = categories[name];

                        }

                        let colors = [];
                        for (let i = 0, n = values.length; i < n; i++) {
                            let rgb = color(colorScale(values[i]));
                            colors.push([rgbScale(rgb.r), rgbScale(rgb.g), rgbScale(rgb.b)]);
                        }

                        // TODO show tooltips, plot higher values on top?
                        let chartData = {
                            hoverinfo: 'text',
                            showlegend: false,
                            name: name,
                            mode: 'markers',
                            type: is3d ? 'scatter3d' : 'scattergl',
                            hoveron: 'points',
                            x: x,
                            y: y,
                            marker: {
                                size: markerSize,
                                color: colors,
                                opacity: markerOpacity,
                                showscale: false,
                            },
                            unselected: {marker: {opacity: unselectedMarkerOpacity, size: unselectedMarkerSize}},
                            selectedpoints: selectedValueCounts.indices && selectedValueCounts.indices.length === 0 ? null : selectedValueCounts.indices,
                            values: values,
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
                    dispatch(_setEmbeddingData(embeddingData.slice(0)));

                },
            );
        }
        Promise.all(promises).finally(() => {
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



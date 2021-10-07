import React from 'react';
import {render} from 'react-dom';
import {Provider} from 'react-redux';
import {applyMiddleware, createStore} from 'redux';
import thunkMiddleware from 'redux-thunk';
import {initGapi, SET_DATASET, SET_EMAIL, SET_SERVER_INFO} from './actions';
import rootReducer from './reducers';
import AppWrapper from './AppWrapper';
// import * as serviceWorker from './serviceWorker';
import mixpanel from 'mixpanel-browser';

let useMixPanel = false;


const logger = store => next => action => {
    if (action.type === SET_SERVER_INFO) {
        if (action.payload.mixpanel) {
            mixpanel.init(action.payload.mixpanel);
            useMixPanel = true;
        }
    }
    if (useMixPanel) {
        if (action.type === SET_DATASET) {
            mixpanel.track('Open Dataset', {name: action.payload.name, id: action.payload.id});
        } else if (action.type === SET_EMAIL) {
            mixpanel.identify(action.payload);
        }
    }
    return next(action);
};
const store = createStore(
    rootReducer,
    applyMiddleware(
        thunkMiddleware, logger
    )
);

function main() {
    render(
        <Provider store={store}>
            {/*<React.StrictMode>*/}

            <AppWrapper/>

            {/*</React.StrictMode>*/}
        </Provider>,
        document.getElementById('root')
    );
}

store.dispatch(initGapi());
main();

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: http://bit.ly/CRA-PWA
// serviceWorker.unregister();

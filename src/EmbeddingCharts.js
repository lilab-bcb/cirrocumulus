import {find} from 'lodash';
import React from 'react';

import {connect} from 'react-redux';
import {getTraceKey, setPrimaryChartSize} from './actions';
import EmbeddingChart from './EmbeddingChart';


class EmbeddingCharts extends React.PureComponent {

    constructor(props) {
        super(props);
        this.resizeListener = () => {
            let width = window.innerWidth - 280;
            let height = Math.max(50, window.innerHeight - 370);
            this.props.handlePrimaryChartSize({width: width, height: height});
            this.windowHeight = window.innerHeight;
        };
        window.addEventListener('resize', this.resizeListener);
    }

    componentWillUnmount() {
        window.removeEventListener('resize', this.resizeListener);
    }


    render() {
        const {activeFeature, embeddingData, onGallery} = this.props;
        if (activeFeature == null) {
            return null;
        }
        const primaryTrace = find(embeddingData, traceInfo => getTraceKey(traceInfo) === activeFeature.embeddingKey);
        if (primaryTrace == null) {
            return null;
        }

        return (<EmbeddingChart
                onGallery={onGallery}
                traceInfo={primaryTrace}
            />
        );
    }

}

const mapStateToProps = state => {
    return {
        activeFeature: state.activeFeature,
        embeddingData: state.embeddingData,
        primaryChartSize: state.primaryChartSize
    };
};
const mapDispatchToProps = dispatch => {
    return {
        handlePrimaryChartSize: value => {
            dispatch(setPrimaryChartSize(value));
        }
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


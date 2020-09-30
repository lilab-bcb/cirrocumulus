import React from 'react';

import {connect} from 'react-redux';
import {getEmbeddingKey, getTraceKey, setPrimaryChartSize} from './actions';
import EmbeddingChart from './EmbeddingChart';


const emptySet = new Set();

class EmbeddingCharts extends React.PureComponent {

    constructor(props) {
        super(props);
        this.resizeListener = () => {
            let width = window.innerWidth - 280;
            // let height = Math.max(1, this.containerElementRef.offsetHeight);
            let height = Math.max(300, window.innerHeight - 220);
            this.props.handlePrimaryChartSize({width: width, height: height});
        };
        window.addEventListener('resize', this.resizeListener);

    }

    componentWillUnmount() {
        window.removeEventListener('resize', this.resizeListener);
    }


    render() {
        const {primaryTraceKey, embeddingData, selection, onGallery, primaryChartSize} = this.props;
        let primaryTraces = embeddingData.filter(traceInfo => getTraceKey(traceInfo) === primaryTraceKey);
        const primaryTrace = primaryTraces.length === 1 ? primaryTraces[0] : null;
        let userPoints = emptySet;
        if (primaryTrace) {
            const embedding = primaryTrace.embedding;
            const fullName = getEmbeddingKey(embedding);
            const chartSelection = selection != null && selection.chart != null ? selection.chart[fullName] : null;
            userPoints = chartSelection ? chartSelection.userPoints : emptySet;
        }

        if (primaryTrace == null) {
            return <div style={{height: primaryChartSize.height}}></div>;
        }

        return (<EmbeddingChart
                onGallery={onGallery}
                traceInfo={primaryTrace}
                selection={userPoints}
            />
        );
    }

}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        primaryTraceKey: state.primaryTraceKey,
        primaryChartSize: state.primaryChartSize,
        selection: state.selection
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


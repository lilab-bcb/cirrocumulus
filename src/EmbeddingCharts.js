import React from 'react';

import {connect} from 'react-redux';
import {Grid} from 'react-virtualized';
import EmbeddingChart from './EmbeddingChart';
import PlotUtil from './PlotUtil';


class EmbeddingCharts extends React.PureComponent {

    getPlots() {
        const embeddingChartSize = this.props.embeddingChartSize;
        const activeTraces = this.props.embeddingData.filter(traceInfo => traceInfo.active);
        const size = PlotUtil.getSingleEmbeddingChartSize();
        let itemSize = Math.floor(size / embeddingChartSize);
        activeTraces.forEach(traceInfo => {

            if (size !== traceInfo.layout.width) {
                traceInfo.layout = Object.assign({}, traceInfo.layout);
                traceInfo.layout.width = itemSize;
                traceInfo.layout.height = itemSize;
            }

            if (traceInfo.data[0].marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].marker.size = 1e-8;
            }
            if (traceInfo.data[0].unselected.marker.size === 0 && traceInfo.data[0].type !== 'scatter3d') {
                traceInfo.data[0].unselected.marker.size = 1e-8;
            }
        });
        if (activeTraces.length === 0) {
            return null;
        }

        let gridWidth = window.screen.availWidth - 280;
        let ncols = Math.min(embeddingChartSize, 2);
        let gridHeight;
        if (embeddingChartSize === 1) {
            gridHeight = itemSize;
        } else if (embeddingChartSize === 2) {
            gridHeight = Math.min(itemSize * 3, window.screen.availHeight - 190);
        } else if (embeddingChartSize === 3) {
            gridHeight = Math.min(itemSize * 2, window.screen.availHeight - 190);
        }


        let elements = [];
        for (let i = 0; i < activeTraces.length; i++) {
            if (i % embeddingChartSize === 0) {
                elements.push([]);
            }
        }

        function cellRenderer({
                                  columnIndex, // Horizontal (column) index of cell
                                  isScrolling, // The Grid is currently being scrolled
                                  isVisible, // This cell is visible within the grid (eg it is not an overscanned cell)
                                  key, // Unique key within array of cells
                                  parent, // Reference to the parent Grid (instance)
                                  rowIndex, // Vertical (row) index of cell
                                  style, // Style object to be applied to cell (to position it);
                                  // This must be passed through to the rendered cell element.
                              }) {


            let index = rowIndex * ncols + columnIndex;
            if (index >= activeTraces.length) {
                return <div style={style} key={key}/>;
            }
            return (
                <EmbeddingChart style={style} traceInfo={activeTraces[index]}
                                key={key}/>
            );
        }

        // can't render more than 6 charts due to webgl context
        return <Grid
            cellRenderer={cellRenderer}
            columnWidth={itemSize + 300}
            columnCount={ncols}
            useDynamicRowHeight={false}
            height={gridHeight}
            overscanColumnCount={0}
            overscanRowCount={0}
            rowHeight={itemSize}
            rowCount={Math.ceil(activeTraces.length / ncols)}
            width={gridWidth}
        />;

    }

    render() {
        return this.getPlots();
    }
}

const mapStateToProps = state => {
    return {
        embeddingData: state.embeddingData,
        embeddingChartSize: state.embeddingChartSize
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(EmbeddingCharts));


import React from 'react';

import {connect} from 'react-redux';
import {Grid} from 'react-virtualized';
import EmbeddingChart from './EmbeddingChart';
import PlotUtil from './PlotUtil';


class EmbeddingCharts extends React.PureComponent {

    getPlots() {
        const embeddingChartSize = this.props.embeddingChartSize;
        const activeTraces = this.props.embeddingData.filter(traceInfo => traceInfo.active);
        const singleChartSize = PlotUtil.getSingleEmbeddingChartSize();
        let itemSize = Math.floor(singleChartSize / embeddingChartSize);
        activeTraces.forEach(traceInfo => {

            if (itemSize !== traceInfo.layout.width) {
                traceInfo.layout = Object.assign({}, traceInfo.layout); // force re-render
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

        // return activeTraces.map(traceInfo => <EmbeddingChart
        //     style={{display: 'inline-block', border: '1px solid LightGrey'}} traceInfo={traceInfo}
        //     key={traceInfo.name + '_' + getEmbeddingKey(traceInfo.data[0].embedding)}/>);

        let gridWidth = window.screen.availWidth - 280;
        let gridHeight = window.screen.availHeight - 190;

        let columnWidth = itemSize + 300; // leave room for legend

        // can't render more than 8 charts due to webgl context
        gridHeight = itemSize * 2;
        let elements = [];
        let row;

        for (let i = 0; i < activeTraces.length; i++) {
            if (i % embeddingChartSize === 0) {
                row = [];
                elements.push(row);
            }
            row.push(activeTraces[i]);
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


            let row = elements[rowIndex];
            let item = row[columnIndex];
            // if (!item || !isVisible || isScrolling) {
            //     return <div style={style} key={key}/>;
            // }
            if (!item) {
                return <div style={style} key={key}/>;
            }
            return (
                <EmbeddingChart style={style} traceInfo={item} key={key}/>
            );
        }


        return <Grid
            cellRenderer={cellRenderer}
            columnWidth={columnWidth}
            columnCount={embeddingChartSize}
            useDynamicRowHeight={false}
            height={gridHeight}
            overscanColumnCount={0}
            overscanRowCount={0}
            rowHeight={itemSize}
            rowCount={elements.length}
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


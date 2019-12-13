import {scaleLinear} from 'd3-scale';
import React from 'react';
import {intFormat, numberFormat} from './formatters';

class CategoricalLegend extends React.PureComponent {


    handleClick = (value, index, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            this.props.handleClick({name: this.props.name, value: value, shiftKey: e.shiftKey, metaKey: e.metaKey});
        }
    };


    render() {
        const {scale, datasetFilter, name, featureSummary, maxHeight, globalFeatureSummary, nObs, nObsSelected} = this.props;
        let clickEnabled = this.props.clickEnabled;
        const categoricalFilterValues = datasetFilter[name];
        const selectionSummary = featureSummary[name];
        let selectedDimensionToCount = {};
        if (selectionSummary != null) {
            for (let i = 0; i < selectionSummary.counts.length; i++) {
                selectedDimensionToCount[selectionSummary.categories[i]] = selectionSummary.counts[i];
            }
        }
        const globalDimensionSummary = globalFeatureSummary[name];
        // if (globalDimensionSummary.max == null) {
        //     let max = 0;
        //     let min = Number.MAX_VALUE;
        //     for (let i = 0; i < globalDimensionSummary.counts.length; i++) {
        //         max = Math.max(max, globalDimensionSummary[i]);
        //         min = Math.min(min, globalDimensionSummary[i]);
        //     }
        //     globalDimensionSummary.min = min;
        //     globalDimensionSummary.max = max;
        // }
        const categories = globalDimensionSummary.categories;
        clickEnabled = clickEnabled && categories.length > 1;

        let maxSize = 60;
        const fractionScale = scaleLinear().domain([0, 1]).range([0, maxSize]).clamp(true);
        return (
            <div className="cirro-chart-legend" style={{
                maxHeight: maxHeight
            }}>
                <b>{name}</b> <small>({categories.length})</small>
                <table>
                    <thead>
                    <tr>
                        {clickEnabled && <td></td>}
                        <td></td>
                        <td><small>{'all'}</small></td>
                        <td><small>{selectionSummary != null ? 'selection' : null}</small></td>
                    </tr>
                    </thead>
                    <tbody>
                    {categories.map((category, i) => {

                        const opacity = categoricalFilterValues == null || categoricalFilterValues.indexOf(category) !== -1 ? 1 : 0.4;
                        // const fractionUnselected = unselectedCounts != null ? unselectedCounts[i] / unselectedTotal : null;
                        // const unselectedSize = unselectedCounts == null ? 0 : fractionScale(fractionUnselected);
                        // const unselectedTitle = unselectedCounts == null ? null : intFormat(unselectedCounts[i]) + ' / ' + intFormat(unselectedTotal) + (unselectedCounts[i] > 0 ? (' ( ' + numberFormat(100 * fractionUnselected) + '%)') : '');
                        const count = globalDimensionSummary.counts[i];
                        const selectedCount = selectedDimensionToCount[category] || 0;
                        const fractionSelected = selectionSummary == null ? 0 : selectedCount / nObsSelected;
                        const selectedSize = fractionScale(fractionSelected);

                        const globalSize = fractionScale(count / nObs);
                        const globalTitle = intFormat(count) + ' / ' + nObs + (' (' + numberFormat(100 * count / nObs) + '%)');
                        const selectionTitle = selectionSummary == null ? null : intFormat(selectedCount) + ' / ' + intFormat(nObsSelected) + (selectedCount > 0 ? (' (' + numberFormat(100 * fractionSelected) + '%)') : '');
                        return <tr
                            style={{cursor: clickEnabled ? 'pointer' : null, opacity: opacity}}
                            onClick={(e) => this.handleClick(category, i, e)} key={category}>
                            {clickEnabled && <td>
                                <div style={{
                                    display: 'inline-block',
                                    width: 10,
                                    height: 10,
                                    background: scale(category)
                                }}/>
                            </td>}
                            <td>
                                <div style={{
                                    maxWidth: 140,
                                    textOverflow: 'ellipsis',
                                    overflow: 'hidden',
                                    display: 'inline-block',
                                    userSelect: 'none'
                                }} title={'' + category}>{'' + category}</div>
                            </td>

                            <td>
                                <div
                                    title={globalTitle}
                                    style={{
                                        display: 'inline-block',
                                        position: 'relative',
                                        width: maxSize,
                                        border: '1px solid black',
                                        height: 9
                                    }}>

                                    <div style={{
                                        position: 'absolute',
                                        width: globalSize,
                                        left: 0,
                                        top: 0,
                                        backgroundColor: 'LightGrey',
                                        height: 9
                                    }}/>
                                </div>
                            </td>
                            {selectionSummary && <td>
                                <div
                                    title={selectionTitle}
                                    style={{
                                        display: 'inline-block',
                                        position: 'relative',
                                        width: maxSize,
                                        border: '1px solid black',
                                        height: 9
                                    }}>

                                    <div style={{
                                        position: 'absolute',
                                        width: selectedSize,
                                        left: 0,
                                        top: 0,
                                        backgroundColor: 'LightGrey',
                                        height: 9
                                    }}/>
                                </div>
                            </td>}
                        </tr>;
                    })
                    }</tbody>

                </table>
            </div>);
    }
}

export default CategoricalLegend;



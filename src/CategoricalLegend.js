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
        const {scale, categoricalFilter, name, featureSummary, maxHeight, clickEnabled} = this.props;
        const categoricalFilterValues = categoricalFilter[name];
        const dimensionSummary = featureSummary.dimensions[name];
        const domain = dimensionSummary.categories;
        let selectedCounts = dimensionSummary.selected_counts;
        let unselectedCounts = dimensionSummary.unselected_counts;
        if (selectedCounts == null && unselectedCounts != null) {
            selectedCounts = unselectedCounts;
        }
        let selectedTotal = 0;
        selectedCounts.forEach(value => selectedTotal += value);

        let unselectedTotal = 0;
        let maxSize = 60;
        const fractionScale = scaleLinear().domain([0, 1]).range([0, maxSize]).clamp(true);
        if (unselectedCounts != null) {
            unselectedCounts.forEach(value => unselectedTotal += value);
        }


        return (
            <div style={{
                display: 'inline-block',
                padding: 10,
                verticalAlign: 'top',
                maxHeight: maxHeight,
                overflow: 'auto'
            }}>
                <b>{name}</b>
                <table>
                    <tbody>
                    {domain.map((d, i) => {

                        const opacity = categoricalFilterValues == null || categoricalFilterValues.indexOf(d) !== -1 ? 1 : 0.4;
                        const fractionUnselected = unselectedCounts != null ? unselectedCounts[i] / unselectedTotal : null;
                        const unselectedSize = unselectedCounts == null ? 0 : fractionScale(fractionUnselected);
                        const unselectedTitle = unselectedCounts == null ? null : intFormat(unselectedCounts[i]) + ' / ' + intFormat(unselectedTotal) + (unselectedCounts[i] > 0 ? (' ( ' + numberFormat(100 * fractionUnselected) + '%)') : '');


                        const fractionSelected = selectedCounts[i] / selectedTotal;
                        const selectedSize = fractionScale(fractionSelected);
                        const selectionTitle = intFormat(selectedCounts[i]) + ' / ' + intFormat(selectedTotal) + (selectedCounts[i] > 0 ? (' ( ' + numberFormat(100 * fractionSelected) + '%)') : '');
                        return <tr
                            style={{cursor: clickEnabled ? 'pointer' : null, opacity: opacity}}
                            onClick={(e) => this.handleClick(d, i, e)} key={d}>
                            {clickEnabled && <td>
                                <div style={{
                                    display: 'inline-block',
                                    width: 10,
                                    height: 10,
                                    background: scale(d)
                                }}/>
                            </td>}
                            <td>
                                <div style={{
                                    maxWidth: 140,
                                    textOverflow: 'ellipsis',
                                    overflow: 'hidden',
                                    display: 'inline-block',
                                    userSelect: 'none'
                                }} title={'' + d}>{'' + d}</div>
                            </td>

                            <td>
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
                            </td>
                            {unselectedCounts != null && <td>
                                <div title={unselectedTitle} style={{
                                    display: 'inline-block',
                                    position: 'relative',
                                    width: maxSize,
                                    border: '1px solid black',
                                    height: 9
                                }}>
                                    <div style={{
                                        position: 'absolute',
                                        width: unselectedSize,
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
                    <tfoot>
                    <tr>
                        {clickEnabled && <td></td>}
                        <td></td>
                        <td><small>{unselectedCounts != null ? 'selection' : null}</small></td>
                        <td><small>{unselectedCounts != null ? 'rest' : null}</small></td>
                    </tr>
                    </tfoot>
                </table>
            </div>);
    }
}

export default CategoricalLegend;



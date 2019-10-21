import React from 'react';
import {getLegendSizeHelper} from './PlotUtil';


class CategoricalLegend extends React.PureComponent {


    handleClick = (value, index, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            this.props.handleClick({name: this.props.name, value: value, shiftKey: e.shiftKey, metaKey: e.metaKey});
        }
    };


    render() {
        const {scale, categoricalFilter, name, selectedValueCounts, maxHeight, clickEnabled} = this.props;
        const categoricalFilterValues = categoricalFilter[name];
        const domain = this.props.domain != null ? this.props.domain : scale.domain();
        const selectedCountMap = selectedValueCounts.categories != null ? selectedValueCounts.categories[name] : null;
        let maxSize = 60;

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
                        let legend = getLegendSizeHelper(selectedCountMap, scale, i, selectedValueCounts.count);
                        let opacity = categoricalFilterValues == null || categoricalFilterValues.indexOf(d) !== -1 ? 1 : 0.4;
                        let groupSize = legend.percentTotal * maxSize;
                        let selectedSize = legend.percentSelected * maxSize;
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
                            {!isNaN(selectedSize) ?
                                <td>
                                    <div
                                        title={legend.selectionTitle}
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
                                </td> : null}
                            <td>
                                <div title={legend.title} style={{
                                    display: 'inline-block',
                                    position: 'relative',
                                    width: maxSize,
                                    border: '1px solid black',
                                    height: 9
                                }}>
                                    <div style={{
                                        position: 'absolute',
                                        width: groupSize,
                                        left: 0,
                                        top: 0,
                                        backgroundColor: 'LightGrey',
                                        height: 9
                                    }}/>
                                </div>
                            </td>
                        </tr>;
                    })
                    }</tbody>
                    <tfoot>
                    <tr>
                        {clickEnabled && <td></td>}
                        <td></td>
                        {selectedCountMap != null ?
                            <td><small>selection</small></td> : null}
                        <td><small>{selectedCountMap != null ? 'group' : null}</small></td>
                    </tr>
                    </tfoot>
                </table>
            </div>);
    }
}

export default CategoricalLegend;



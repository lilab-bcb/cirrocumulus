import React from 'react';
import {intFormat, numberFormat} from './formatters';


class ContinuousLegend extends React.PureComponent {


    render() {
        const {name, selectedValueCounts, maxHeight, summary, nTotal} = this.props;
        const selectionSummary = selectedValueCounts.summary != null ? selectedValueCounts.summary[name] : null;
        let summaries = [];
        let totals = [nTotal];
        let names = [''];
        if (summary) {
            summaries.push(summary);
        }
        if (selectionSummary) {
            names[0] = 'All';
            names.push('Selection');
            summaries.push(selectionSummary);
            totals.push(selectedValueCounts.count);
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
                {/*{summaries.length > 0 && <Plot*/}
                {/*    data={chartData.data}*/}
                {/*    layout={chartData.layout}*/}
                {/*    config={chartData.config}*/}
                {/*/>}*/}
                {summaries.length > 0 && <table>
                    <thead>
                    <tr>
                        <td></td>
                        <td>Mean</td>
                        <td>Variance</td>
                        <td>% Expressed</td>
                    </tr>
                    </thead>
                    <tbody>
                    {summaries.map((summary, i) => {
                        return <tr key={i}>
                            <td>{names[i]}</td>
                            <td>{numberFormat(summary.mean)}</td>
                            <td>{numberFormat(summary.variance)}</td>
                            <td>{intFormat(Math.round(100 * (summary.num_expressed / totals[i])))}</td>
                        </tr>;
                    })}
                    </tbody>
                </table>}
            </div>);
    }
}

export default ContinuousLegend;



import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import React from 'react';
import {intFormat, numberFormat} from './formatters';

class ContinuousLegend extends React.PureComponent {


    getChartData(summaries) {
        let layout = {
            xaxis: {
                showbackground: false,
                autorange: true,
                showgrid: false,
                zeroline: false,
                showline: false,
                title: '',

            },
            yaxis: {
                showbackground: false,
                autorange: true,
                showgrid: false,
                zeroline: false,
                showline: false,
                title: '',

            },
            dragmode: false,
            hovermode: 'x',
            width: 180,
            height: 70,
            margin: {
                l: summaries.length > 1 ? 70 : 0,
                b: 0,
                r: 0,
                t: 0,
                autoexpand: true
            },

            autosize: true,
            displaylogo: false,
            showlegend: false,
        };
        //  let sizeScale = scaleLinear().domain([sizeMin, sizeMax]).range([1, maxDiameter]).clamp(true);
        let y;
        if (summaries.length > 1) {
            y = ['All', 'Selection'];
        } else if (summaries.length === 1) {
            y = [''];
        }


        let data = [
            {
                hoverinfo: 'none',
                x: summaries.map(summary => summary.mean),
                y: y,
                sizemode: 'diameter',
                marker: {shape: 'circle'},
                // marker: {size: [10, 20]},
                error_x: {
                    type: 'data',
                    array: summaries.map(summary => summary.variance)
                },
                type: 'scatter',
                mode: 'markers'
            }
        ];
        return {data: data, layout: layout, config: {displayModeBar: false}};
    }

    render() {
        const {name, selectedValueCounts, maxHeight, summary, nTotal} = this.props;
        const selectionSummary = selectedValueCounts.summary != null ? selectedValueCounts.summary[name] : null;
        let summaries = [];
        let totals = [nTotal];
        if (summary) {
            summaries.push(summary);
        }
        if (selectionSummary) {
            summaries.push(selectionSummary);
            totals.push(selectedValueCounts.count);
        }

        let chartData = this.getChartData(summaries);

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
                                <td>{chartData.data[0].y[i]}</td>
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



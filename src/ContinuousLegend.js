import React from 'react';
import {intFormat, numberFormat} from './formatters';


class ContinuousLegend extends React.PureComponent {


    getValues(measureSummary, statistic) {
        if (statistic === 'Mean') {
            return measureSummary.mean.map((value, i) => {
                return <td key={i}>{numberFormat(value)}</td>;
            });

        } else if (statistic === '% Expressed') {
            return measureSummary.fraction_expressed.map((value, i) => {
                return <td key={i}>{intFormat(100 * value)}</td>;
            });
        }
    }

    render() {
        const {name, featureSummary, maxHeight} = this.props;
        const displayName = name === '__count' ? 'count' : name;
        const measureSummary = featureSummary.measures[name];
        const summaryNames = measureSummary != null && measureSummary.mean.length === 2 ? ['selection', 'rest'] : [''];
        const statistics = ['Mean', '% Expressed'];

        return (
            <div style={{
                display: 'inline-block',
                padding: 10,
                verticalAlign: 'top',
                maxHeight: maxHeight,
                overflow: 'auto'
            }}>
                <b>{displayName}</b>

                {measureSummary && <table>

                    <tbody>
                    {statistics.map((statistic, i) => {
                        return <tr key={i}>
                            <td style={{textAlign: 'right'}}>{statistic}:</td>
                            {this.getValues(measureSummary, statistic)}
                        </tr>;
                    })}
                    </tbody>
                    <tfoot>
                    <tr>
                        <td></td>
                        <td><small>{summaryNames.length === 2 ? 'selection' : null}</small></td>
                        <td><small>{summaryNames.length === 2 ? 'rest' : null}</small></td>
                    </tr>
                    </tfoot>
                </table>}
            </div>);
    }
}

export default ContinuousLegend;



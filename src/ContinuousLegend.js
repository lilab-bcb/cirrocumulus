import React from 'react';
import {numberFormat0, numberFormat} from './formatters';


class ContinuousLegend extends React.PureComponent {


    getSelectedAndUnselectedTds(measureSummary, statistic) {
        if (statistic === 'Mean') {
            return measureSummary.mean.map((value, i) => {
                return <td key={i}>{numberFormat0(value)}</td>;
            });

        } else if (statistic === '% Expressed') {
            return measureSummary.fraction_expressed.map((value, i) => {
                return <td key={i}>{numberFormat0(100 * value)}</td>;
            });
        }
    }

    getGlobalSummaryTd(measureSummary, statistic) {
        if (statistic === 'Mean') {
            return <td>{numberFormat(measureSummary.mean)}</td>;
        } else if (statistic === '% Expressed') {
            return <td>{}</td>;
        }
    }

    getTable(statistics, summaryNames, selectionSummary, globalSummary) {
        return (<table>
            <tbody>
            <tr>
                <td style={{textAlign: 'right'}}>{'Mean'}:</td>
                <td>{numberFormat(globalSummary.mean)}</td>
                {selectionSummary && <td>{numberFormat(selectionSummary.mean)}</td>}
            </tr>
            <tr>
                <td style={{textAlign: 'right'}}>{'% Expressed'}:</td>
                <td>{numberFormat0(100 * globalSummary.num_expressed / this.props.nObs)}</td>
                {selectionSummary &&
                <td>{numberFormat0(100 * selectionSummary.num_expressed / this.props.nObsSelected)}</td>}
            </tr>
            </tbody>
            <tfoot>
            <tr>
                <td></td>
                {summaryNames.map(summaryName => {
                    return <td key={summaryName}><small>{summaryName}</small></td>;
                })}
            </tr>
            </tfoot>
        </table>);
    }

    render() {
        const {name, featureSummary, globalFeatureSummary, maxHeight} = this.props;
        const displayName = name === '__count' ? 'count' : name;
        const selectionSummary = featureSummary[name];
        const globalSummary = globalFeatureSummary[name];
        const summaryNames = ['all'];
        // TODO compute unselected mean and % expressed from globals
        if (selectionSummary != null) {
            summaryNames.push('selection');
            //  summaryNames.push('rest');
        }
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
                {globalSummary != null && name !== '__count' && this.getTable(statistics, summaryNames, selectionSummary, globalSummary)}
            </div>);
    }
}

export default ContinuousLegend;



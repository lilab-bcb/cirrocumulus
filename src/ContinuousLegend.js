import React from 'react';
import {numberFormat, numberFormat0} from './formatters';


class ContinuousLegend extends React.PureComponent {


    getTable(statistics, summaryNames, selectionSummary, globalSummary) {
        return (<table>
            <thead>
            <tr>
                <td></td>
                {summaryNames.map(summaryName => {
                    return <td key={summaryName}><small>{summaryName}</small></td>;
                })}
            </tr>
            </thead>
            <tbody>
            <tr>
                <td style={{textAlign: 'right'}}>{'Mean'}:</td>
                <td>{numberFormat(globalSummary.mean)}</td>
                {selectionSummary && <td>{numberFormat(selectionSummary.mean)}</td>}
            </tr>
            {globalSummary.numExpressed != null && <tr>
                <td style={{textAlign: 'right'}}>{'% Expressed'}:</td>
                <td>{numberFormat0(100 * globalSummary.numExpressed / this.props.nObs)}</td>
                {selectionSummary &&
                <td>{numberFormat0(100 * selectionSummary.numExpressed / this.props.nObsSelected)}</td>}
            </tr>}
            </tbody>

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



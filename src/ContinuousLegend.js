import React from 'react';
import {numberFormat, numberFormat0} from './formatters';
import {stripTrailingZeros} from './util';


class ContinuousLegend extends React.PureComponent {


    getTable(summaryNames, selectionSummary, globalSummary) {
        return (
            <table>
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
                    <td>{stripTrailingZeros(numberFormat(globalSummary.mean))}</td>
                    {selectionSummary && <td>{stripTrailingZeros(numberFormat(selectionSummary.mean))}</td>}
                </tr>

                {globalSummary.numExpressed != null &&
                <tr>
                    <td style={{textAlign: 'right'}}>{'% Expressed'}:</td>
                    <td>{numberFormat0(100 * globalSummary.numExpressed / this.props.nObs)}</td>
                    {selectionSummary &&
                    <td>{numberFormat0(100 * selectionSummary.numExpressed / this.props.nObsSelected)}</td>}
                </tr>}
                {(globalSummary.numExpressed == null || globalSummary.min !== 0) && <tr>
                    <td style={{textAlign: 'right'}}>{'Min'}:</td>
                    <td>{stripTrailingZeros(numberFormat(globalSummary.min))}</td>
                    {selectionSummary && <td>{stripTrailingZeros(numberFormat(selectionSummary.min))}</td>}
                </tr>}
                <tr>
                    <td style={{textAlign: 'right'}}>{'Max'}:</td>
                    <td>{stripTrailingZeros(numberFormat(globalSummary.max))}</td>
                    {selectionSummary && <td>{stripTrailingZeros(numberFormat(selectionSummary.max))}</td>}
                </tr>
                </tbody>
            </table>
        );
    }

    render() {
        const {name, featureSummary, globalFeatureSummary} = this.props;
        const selectionSummary = featureSummary[name];
        const globalSummary = globalFeatureSummary[name];
        const summaryNames = ['all'];
        // TODO compute unselected mean and % expressed from globals
        if (selectionSummary != null) {
            summaryNames.push('selection');
            //  summaryNames.push('rest');
        }

        return (
            <>
                {globalSummary != null && name !== '__count' && this.getTable(summaryNames, selectionSummary, globalSummary)}
            </>);
    }
}

export default ContinuousLegend;



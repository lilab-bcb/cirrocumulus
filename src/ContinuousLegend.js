import React from 'react';
import {numberFormat, numberFormat0} from './formatters';
import {FEATURE_TYPE, stripTrailingZeros} from './util';


function ContinuousLegend(props) {

    const {name, featureSummary, globalFeatureSummary, nObs, nObsSelected, type} = props;
    const selectionSummary = featureSummary[name];
    const globalSummary = globalFeatureSummary[name];
    const summaryNames = ['all'];
    // TODO compute unselected mean and % expressed from globals
    if (selectionSummary != null) {
        summaryNames.push('selection');
        //  summaryNames.push('rest');
    }

    function getTable(summaryNames, selectionSummary, globalSummary) {
        const showPercentExpressed = globalSummary.numExpressed != null && type === FEATURE_TYPE.X;
        const showMin = type !== FEATURE_TYPE.X;
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
                {showPercentExpressed &&
                <tr>
                    <td style={{textAlign: 'right'}}>{'% Expressed'}:</td>
                    <td>{numberFormat0(100 * globalSummary.numExpressed / nObs)}</td>
                    {selectionSummary &&
                    <td>{numberFormat0(100 * selectionSummary.numExpressed / nObsSelected)}</td>}
                </tr>}
                {showMin && <tr>
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


    return (
        <>
            {globalSummary != null && name !== '__count' && getTable(summaryNames, selectionSummary, globalSummary)}
        </>);

}

export default ContinuousLegend;



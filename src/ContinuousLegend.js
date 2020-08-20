import React from 'react';
import {numberFormat, numberFormat0} from './formatters';


class ContinuousLegend extends React.PureComponent {


    getTable(summaryNames, selectionSummary, globalSummary) {
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
            {(globalSummary.numExpressed == null || globalSummary.min !== 0) && <tr>
                <td style={{textAlign: 'right'}}>{'Min'}:</td>
                <td>{numberFormat(globalSummary.min)}</td>
                {selectionSummary && <td>{numberFormat(selectionSummary.min)}</td>}
            </tr>}
            <tr>
                <td style={{textAlign: 'right'}}>{'Max'}:</td>
                <td>{numberFormat(globalSummary.max)}</td>
                {selectionSummary && <td>{numberFormat(selectionSummary.max)}</td>}
            </tr>

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
        return (
            <React.Fragment>
                <h4 style={{margin: 0}}>{displayName}</h4>
                {globalSummary != null && name !== '__count' && this.getTable(summaryNames, selectionSummary, globalSummary)}
            </React.Fragment>);
    }
}

export default ContinuousLegend;



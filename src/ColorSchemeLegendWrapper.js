import React from 'react';
import ContinuousLegend from './ContinuousLegend';
import MeasureFilter from './MeasureFilter';

export default function ColorSchemeLegendWrapper(props) {
    const {
        datasetFilter,
        featureSummary,
        globalFeatureSummary,
        handleUpdate,
        name,
        nObs,
        nObsSelected,
        style,
        type
    } = props;

    const legendStyle = Object.assign({display: 'inline-block', verticalAlign: 'top'}, style || {});
    const isCount = name === '__count';
    return (
        <div data-testid="continuous-legend" className="cirro-condensed" style={legendStyle}>
            {/*ContinuousLegend shows stats table */}
            <ContinuousLegend name={name}
                              featureSummary={featureSummary}
                              nObs={nObs}
                              nObsSelected={nObsSelected}
                              globalFeatureSummary={globalFeatureSummary}
                              type={type}></ContinuousLegend>
            {/*{!isCount && handleDomain &&*/}
            {/*<InputLabel shrink={true} variant={"standard"}>Custom Color Scale</InputLabel>}*/}
            {/*{name !== '__count' && handleDomain &&*/}
            {/*<TextField InputLabelProps={{shrink: true}} margin="none"*/}
            {/*           style={{maxWidth: 60, marginRight: 4, marginTop: 0}}*/}
            {/*           size="small" type="text"*/}

            {/*           onChange={onMinChange} label="Min"*/}
            {/*           value={min}/>}*/}
            {/*{!isCount && handleDomain &&*/}
            {/*<TextField InputLabelProps={{shrink: true}} margin="none" style={{maxWidth: 60, marginTop: 0}}*/}
            {/*           size="small"*/}
            {/*           type="text"*/}

            {/*           onChange={onMaxChange} label="Max"*/}
            {/*           value={max}/>}*/}
            {!isCount &&
            <MeasureFilter datasetFilter={datasetFilter} name={name} handleUpdate={handleUpdate}/>}
        </div>);

}


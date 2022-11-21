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
    type,
  } = props;

  const legendStyle = Object.assign(
    {display: 'inline-block', verticalAlign: 'top'},
    style || {},
  );
  const isCount = name === '__count';
  return (
    <div
      data-testid="continuous-legend"
      className="cirro-condensed"
      style={legendStyle}
    >
      {/*ContinuousLegend shows stats table */}
      <ContinuousLegend
        name={name}
        featureSummary={featureSummary}
        nObs={nObs}
        nObsSelected={nObsSelected}
        globalFeatureSummary={globalFeatureSummary}
        type={type}
      ></ContinuousLegend>
      {!isCount && (
        <MeasureFilter
          datasetFilter={datasetFilter}
          name={name}
          handleUpdate={handleUpdate}
        />
      )}
    </div>
  );
}

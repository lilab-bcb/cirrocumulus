import React from 'react';

import ColorSchemeLegend from './ColorSchemeLegend';
import ContinuousLegend from './ContinuousLegend';
import MeasureFilter from './MeasureFilter';

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {
        const {scale, name, nObs, nObsSelected, featureSummary, globalFeatureSummary, maxHeight, style, datasetFilter, handleUpdate} = this.props;
        return (
            <div className="cirro-chart-legend" style={{
                maxHeight: maxHeight
            }}>
                <ContinuousLegend name={name}
                                  summary={scale.summary}
                                  featureSummary={featureSummary}
                                  nObs={nObs}
                                  nObsSelected={nObsSelected}
                                  globalFeatureSummary={globalFeatureSummary}
                                  maxHeight={maxHeight}></ContinuousLegend>
                {name !== '__count' &&
                <MeasureFilter datasetFilter={datasetFilter} name={name} handleUpdate={handleUpdate}/>}
                <ColorSchemeLegend width={this.props.width} height={this.props.height} style={style} scale={scale}
                                   label={this.props.label}></ColorSchemeLegend></div>);

    }
}

export default ColorSchemeLegendWrapper;


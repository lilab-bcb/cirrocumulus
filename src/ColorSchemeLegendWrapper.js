import React from 'react';

import ColorSchemeLegend from './ColorSchemeLegend';
import ContinuousLegend from './ContinuousLegend';
import MeasureFilter from './MeasureFilter';

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {
        const {scale, name, nObs, nObsSelected, featureSummary, globalFeatureSummary, maxHeight, datasetFilter, handleUpdate, showColorScheme} = this.props;
        return (
            <div className="cirro-condensed" style={{display: 'inline-block', verticalAlign: 'top'}}>
                {/*ContinuousLegend shows stats table */}
                <ContinuousLegend name={name}
                                  summary={scale.summary}
                                  featureSummary={featureSummary}
                                  nObs={nObs}
                                  nObsSelected={nObsSelected}
                                  globalFeatureSummary={globalFeatureSummary}
                                  maxHeight={maxHeight}></ContinuousLegend>
                {name !== '__count' &&
                <MeasureFilter datasetFilter={datasetFilter} name={name} handleUpdate={handleUpdate}/>}
                {/*{showColorScheme &&*/}
                {/*<ColorSchemeLegend width={this.props.width} height={this.props.height} style={{paddingLeft: 10}}*/}
                {/*                   scale={scale}*/}
                {/*                   label={this.props.label}></ColorSchemeLegend>}*/}
            </div>);

    }
}

export default ColorSchemeLegendWrapper;

import React from 'react';

import ColorSchemeLegend from './ColorSchemeLegend';
import ContinuousLegend from './ContinuousLegend';
import MeasureFilter from './MeasureFilter';

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {
        const {scale, name, featureSummary, maxHeight, nTotal, style, datasetFilter, handleUpdate} = this.props;
        return (
            <div style={{
                display: 'inline-block',
                padding: 10,
                verticalAlign: 'top',
                maxHeight: maxHeight,
                overflow: 'auto'
            }}>
                <ContinuousLegend name={name}
                                  nTotal={nTotal}
                                  summary={scale.summary}
                                  featureSummary={featureSummary}
                                  maxHeight={maxHeight}></ContinuousLegend>
                {name !== '__count' &&
                <MeasureFilter datasetFilter={datasetFilter} name={name} handleUpdate={handleUpdate}/>}
                <ColorSchemeLegend width={this.props.width} height={this.props.height} style={style} scale={scale}
                                   label={this.props.label}></ColorSchemeLegend></div>);

    }
}

export default ColorSchemeLegendWrapper;


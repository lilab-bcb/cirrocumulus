import React from 'react';

import ColorSchemeLegend from './ColorSchemeLegend';
import ContinuousLegend from './ContinuousLegend';

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {
        const {scale, name, featureSummary, maxHeight, nTotal, style} = this.props;
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
                <ColorSchemeLegend width={this.props.width} height={this.props.height} style={style} scale={scale}
                                   label={this.props.label}></ColorSchemeLegend></div>);

    }
}

export default ColorSchemeLegendWrapper;


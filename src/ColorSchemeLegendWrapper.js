import React from 'react';
import CategoricalLegend from './CategoricalLegend';

import ColorSchemeLegend from './ColorSchemeLegend';

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {
        const {scale, name, selectedValueCounts, maxHeight, style} = this.props;
        const selectedCountMap = selectedValueCounts.categories != null ? selectedValueCounts.categories[name] : null;
        const domain = scale.valueCounts.values;
        if (selectedCountMap) {
            // TODO fix hack below
            if (selectedCountMap['True'] != null) {
                selectedCountMap['true'] = selectedCountMap['True'];
            }
            if (selectedCountMap['False'] != null) {
                selectedCountMap['false'] = selectedCountMap['False'];
            }
        }

        return (
            <div style={{
                display: 'inline-block',
                padding: 10,
                verticalAlign: 'top',
                maxHeight: maxHeight,
                overflow: 'auto'
            }}>
                <CategoricalLegend scale={scale} name={name}
                                   selectedValueCounts={selectedValueCounts}
                                   maxHeight={maxHeight}
                                   domain={domain} clickEnabled={false}
                                   legendVisibility={{}}></CategoricalLegend>
                <ColorSchemeLegend width={this.props.width} height={this.props.height} style={style} scale={scale}
                                   label={this.props.label}></ColorSchemeLegend></div>);

    }
}

export default ColorSchemeLegendWrapper;


import {format} from 'd3-format';
import React from 'react';

import ColorSchemeLegend from './ColorSchemeLegend';
import {getLegendSizeHelper, getLegendSizeScale} from './PlotUtil';

const intFormat = format(',');

class ColorSchemeLegendWrapper extends React.PureComponent {

    render() {

        const style = this.props.style;
        const scale = this.props.scale;
        const name = this.props.name;


        let values = scale.valueCounts.values;
        let domain = scale.domain();
        if (domain[0] === domain[1]) {
            values = [];
        }

        const selectedValueCounts = this.props.selectedValueCounts;
        const selectedCountMap = selectedValueCounts.categories != null ? selectedValueCounts.categories[name] : null;
        if (selectedCountMap) {
            // TODO fix hack below
            if (selectedCountMap['True'] != null) {
                selectedCountMap['true'] = selectedCountMap['True'];
            }
            if (selectedCountMap['False'] != null) {
                selectedCountMap['false'] = selectedCountMap['False'];
            }
        }

        let sizeScale = getLegendSizeScale(selectedCountMap, scale.valueCounts.values, scale.valueCounts.counts);


        return (
            <div style={style}>{values.map((d, i) => {
                let label = d;
                if (label == true) {
                    label = 'expressed';
                } else if (label == false) {
                    label = 'not expressed';
                }
                let legend = getLegendSizeHelper(selectedCountMap, scale, sizeScale, i);
                return <div key={d}>

                    <div style={{
                        display: 'inline-block',
                        width: legend.width,
                        height: 10,
                        background: 'LightGrey'
                    }}/>
                    <label
                        style={{marginLeft: 4}}>{label + ' - ' + legend.text}</label>
                </div>;
            })
            }<ColorSchemeLegend width={this.props.width} height={this.props.height} style={style} scale={scale}
                                label={this.props.label}></ColorSchemeLegend></div>);

    }
}

export default ColorSchemeLegendWrapper;


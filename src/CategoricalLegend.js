import React from 'react';
import {getLegendSizeHelper, getLegendSizeScale} from './PlotUtil';


class CategoricalLegend extends React.PureComponent {

    constructor(props) {
        super(props);
    }

    render() {
        const scale = this.props.scale;
        const domain = scale.domain();

        const selectedValueCounts = this.props.selectedValueCounts;
        const selectedCountMap = selectedValueCounts.categories != null ? selectedValueCounts.categories[this.props.name] : null;
        let sizeScale = getLegendSizeScale(selectedCountMap, scale.valueCounts.values, scale.valueCounts.counts);
        return (
            <div style={{display: 'inline-block', padding: 10, verticalAlign: 'top'}}>{domain.map((d, i) => {
                let legend = getLegendSizeHelper(selectedCountMap, scale, sizeScale, i);
                return <div key={d}>

                    <div style={{
                        display: 'inline-block',
                        width: legend.width,
                        height: 10,
                        background: scale(d)
                    }}/>
                    <label
                        style={{marginLeft: 4}}>{d + ' - ' + legend.text}</label>
                </div>;
            })
            }</div>);
    }
}

export default CategoricalLegend;



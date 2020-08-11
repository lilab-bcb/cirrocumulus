import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotSortOrder} from './actions';
import DotPlotCanvas from './DotPlotCanvas';


class DotPlots extends React.PureComponent {
    render() {
        const {dotPlotData, categoricalNames, selectedDotPlotData} = this.props;

        if (dotPlotData.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        let selectedDotPlotNameToData = {};
        selectedDotPlotData.forEach(categoryItem => {
            selectedDotPlotNameToData[categoryItem.name] = categoryItem;
        });
        return <div>
            {dotPlotData.map((categoryItem) => {
                let selectedData = selectedDotPlotNameToData[categoryItem.name];
                let meanRange = categoryItem.meanRange;
                let fractionRange = categoryItem.fractionRange;
                if (selectedData != null) {
                    meanRange = meanRange.slice();
                    fractionRange = fractionRange.slice();
                    meanRange[0] = Math.min(meanRange[0], selectedData.meanRange[0]);
                    meanRange[1] = Math.max(meanRange[1], selectedData.meanRange[1]);
                    fractionRange[0] = Math.min(fractionRange[0], selectedData.fractionRange[0]);
                    fractionRange[1] = Math.max(fractionRange[1], selectedData.fractionRange[1]);

                }
                if (meanRange[0] === meanRange[1]) {
                    meanRange[1]++;
                }
                if (fractionRange[0] === fractionRange[1]) {
                    fractionRange[0] = 0;
                    fractionRange[1] = 1;
                }
                return (
                    <React.Fragment key={categoryItem.name}>
                        <DotPlotCanvas renamedCategories={categoricalNames[categoryItem.name]}
                                       legend={selectedData == null}
                                       meanRange={meanRange}
                                       fractionRange={fractionRange}
                                       onSortOrderChanged={this.props.onSortOrderChanged}
                                       data={categoryItem}/>
                        {selectedData != null ?
                            <DotPlotCanvas renamedCategories={categoricalNames[categoryItem.name]}
                                           subtitle="selection"
                                           legend={true}
                                           meanRange={meanRange}
                                           fractionRange={fractionRange}
                                           data={selectedData}/> : null}
                    </React.Fragment>
                );

            })}
        </div>;
    }
}

const mapStateToProps = state => {
    return {
        dotPlotData: state.dotPlotData,
        selectedDotPlotData: state.selectedDotPlotData,
        categoricalNames: state.categoricalNames
    };
};
const mapDispatchToProps = dispatch => {
    return {
        onSortOrderChanged: (payload) => {
            dispatch(setDotPlotSortOrder(payload));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DotPlots));


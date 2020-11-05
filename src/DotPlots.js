import React from 'react';

import {connect} from 'react-redux';
import {setDotPlotInterpolator, setDotPlotSortOrder} from './actions';
import {DotPlotGroup} from './DotPlotGroup';

class DotPlots extends React.PureComponent {
    render() {
        const {chartOptions, dotPlotData, dotPlotInterpolator, categoricalNames, handleInterpolator, onSortOrderChanged, selectedDotPlotData} = this.props;

        if (dotPlotData.length === 0) {
            return <h4>Please enter one or more categorical observations and one or more features.</h4>;
        }
        const textColor = chartOptions.darkMode ? 'white' : 'black';
        let selectedDotPlotNameToData = {};
        selectedDotPlotData.forEach(categoryItem => {
            selectedDotPlotNameToData[categoryItem.name] = categoryItem;
        });
        return <div>
            {dotPlotData.map((categoryItem) => {
                let selectedData = selectedDotPlotNameToData[categoryItem.name];
                let renamedCategories = categoricalNames[categoryItem.name];
                return <DotPlotGroup key={categoryItem.name}
                                     categoryItem={categoryItem}
                                     selectedData={selectedData}
                                     interpolator={dotPlotInterpolator}
                                     handleInterpolator={handleInterpolator}
                                     onSortOrderChanged={onSortOrderChanged}
                                     renamedCategories={renamedCategories}
                                     textColor={textColor}/>;
            })}
        </div>;
    }
}

const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions,
        dotPlotInterpolator: state.dotPlotInterpolator,
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
        handleInterpolator: value => {
            dispatch(setDotPlotInterpolator(value));
        },
    };
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(DotPlots));


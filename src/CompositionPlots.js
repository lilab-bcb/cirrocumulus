import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import React from 'react';
import {connect} from 'react-redux';
import CompositionPlot from './CompositionPlot';
import {FEATURE_TYPE, NATSORT} from './util';

function getColorScale(embeddingData, dimension) {
    let categoryColorScale = null;
    for (let i = 0; i < embeddingData.length; i++) {
        if (dimension === embeddingData[i].name) {
            categoryColorScale = embeddingData[i].colorScale; // TODO make color scale independent of embedding
            break;
        }
    }
    if (categoryColorScale == null) {
        categoryColorScale = scaleOrdinal(schemeCategory10); // TODO make color scale independent of embedding
    }
    return categoryColorScale;
}

function createSorter(name, categoryOrder, categoricalNames) {
    if (categoryOrder[name]) {
        const orderedCategories = categoryOrder[name];
        const categoryToIndex = new Map();
        const nameMap = categoricalNames[name] || {};
        for (let i = 0; i < orderedCategories.length; i++) {
            let value = orderedCategories[i];
            let newValue = nameMap[value];
            if (newValue !== undefined) {
                value = newValue;
            }
            categoryToIndex.set(value, i);
        }
        return (a, b) => categoryToIndex.get(a) - categoryToIndex.get(b);
    }
    return NATSORT;

}

function getComposition(dataset, obsCat, cachedData, categoricalNames, selection) {
    const ncategories = obsCat.length;

    if (ncategories >= 2) {
        let categoryValues = [];
        let nObs = dataset.shape[0];
        const renamedDimensions = [];
        for (let categoryIndex = 0; categoryIndex < ncategories; categoryIndex++) {
            const array = cachedData[obsCat[categoryIndex]];
            if (array == null) {
                return null;
            }
            categoryValues.push(array);
            renamedDimensions.push(categoricalNames[obsCat[categoryIndex]] || {});
        }
        const hasSelection = selection != null && selection.size > 0;
        const dimensionIndex = ncategories - 1;
        const seriesToValueToCounts = {};
        const seriesToSeriesArray = {};
        for (let i = 0; i < nObs; i++) {
            if (hasSelection && !selection.has(i)) {
                continue;
            }
            const seriesArray = [];
            for (let categoryIndex = 0; categoryIndex < ncategories - 1; categoryIndex++) {
                let value = categoryValues[categoryIndex][i];
                const nameMap = renamedDimensions[categoryIndex];
                let newValue = nameMap[value];
                if (newValue !== undefined) {
                    value = newValue;
                }
                seriesArray.push(value);
            }
            const seriesKey = seriesArray.join(',');
            let valueToCounts = seriesToValueToCounts[seriesKey];
            if (valueToCounts === undefined) {
                valueToCounts = {};
                seriesToSeriesArray[seriesKey] = seriesArray;
                seriesToValueToCounts[seriesKey] = valueToCounts;
            }
            let category = categoryValues[dimensionIndex][i];
            const nameMap = renamedDimensions[dimensionIndex];
            let newCategory = nameMap[category];
            if (newCategory !== undefined) {
                category = newCategory;
            }
            const count = valueToCounts[category] || 0;
            valueToCounts[category] = count + 1;
        }

        const sorters = [];
        const categoryOrder = dataset.categoryOrder || {};
        for (let categoryIndex = 0; categoryIndex < ncategories - 1; categoryIndex++) {
            const name = obsCat[categoryIndex];
            sorters.push(createSorter(name, categoryOrder, categoricalNames));
        }
        const series = Object.keys(seriesToValueToCounts);

        series.sort((a, b) => {
            const val1 = seriesToSeriesArray[a];
            const val2 = seriesToSeriesArray[b];
            for (let i = 0; i < sorters.length; i++) {
                const r = sorters[i](val1[i], val2[i]);
                if (r !== 0) {
                    return r;
                }
            }
        });

        let uniqueValuesSet = new Set();
        for (let key in seriesToValueToCounts) {
            const valueToCounts = seriesToValueToCounts[key];
            for (const value in valueToCounts) {
                uniqueValuesSet.add(value);
            }
        }

        const uniqueValues = Array.from(uniqueValuesSet);
        uniqueValues.sort(createSorter(obsCat[ncategories - 1], categoryOrder, categoricalNames));
        return {seriesToValueToCounts: seriesToValueToCounts, uniqueValues: uniqueValues, series: series};
    }
    return null;
}


function CompositionPlots(props) {
    const {cachedData, categoricalNames, chartOptions, dataset, embeddingData, searchTokens, selection} = props;
    const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT).map(item => item.value);
    if (obsCat.length > 1) {
        const dimension = obsCat[obsCat.length - 1];
        const colorScale = getColorScale(embeddingData, dimension);
        const composition = getComposition(dataset, obsCat, cachedData, categoricalNames);
        if (composition == null) {
            return null;
        }
        const textColor = chartOptions.darkMode ? 'white' : 'black';
        const selectedComposition = selection.size > 0 ? getComposition(dataset, obsCat, cachedData, categoricalNames, selection) : null;
        const title = dimension + ' composition in ' + obsCat.slice(0, obsCat.length - 1).join(', ');
        return <><CompositionPlot seriesToValueToCounts={composition.seriesToValueToCounts}
                                  dimension={dimension}
                                  title={title}
                                  colorScale={colorScale} series={composition.series}
                                  uniqueValues={composition.uniqueValues}
                                  textColor={textColor}/>

            {selectedComposition &&
            <CompositionPlot seriesToValueToCounts={selectedComposition.seriesToValueToCounts}
                             dimension={dimension}
                             title={title}
                             subtitle="selection"
                             colorScale={colorScale} series={selectedComposition.series}
                             uniqueValues={selectedComposition.uniqueValues}
                             textColor={textColor}/>}</>;
    }
    return null;
}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        chartOptions: state.chartOptions,
        dataset: state.dataset,
        embeddingData: state.embeddingData,
        searchTokens: state.searchTokens,
        selection: state.selection
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps,
)(CompositionPlots));


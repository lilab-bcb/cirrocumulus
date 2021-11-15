import {scaleOrdinal} from 'd3-scale';
import {schemeCategory10} from 'd3-scale-chromatic';
import React, {useEffect, useRef, useState} from 'react';
import {connect} from 'react-redux';
import CompositionPlot from './CompositionPlot';
import {FEATURE_TYPE, getCategoryValue, NATSORT} from './util';

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
            value = getCategoryValue(nameMap, value);
            categoryToIndex.set(value, i);
        }
        return (a, b) => categoryToIndex.get(a) - categoryToIndex.get(b);
    }
    return NATSORT;

}

function getComposition(dataset, obsCat, cachedData, categoricalNames, selection) {
    const ncategories = obsCat.length;
    const maxCategories = 100;
    if (ncategories >= 2) {
        let nseries = 0;
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
                value = getCategoryValue(nameMap, value);
                seriesArray.push(value);
            }
            const seriesKey = seriesArray.join(',');
            let valueToCounts = seriesToValueToCounts[seriesKey];
            if (valueToCounts === undefined) {
                valueToCounts = {};
                seriesToSeriesArray[seriesKey] = seriesArray;
                seriesToValueToCounts[seriesKey] = valueToCounts;
                nseries++;
                if (nseries > maxCategories) {
                    return null;
                }
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
            return 0;
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
    const {
        cachedData,
        categoricalNames,
        darkMode,
        dataset,
        embeddingData,
        searchTokens,
        selection,
        visible
    } = props;
    const [needsUpdate, setNeedsUpdate] = useState(true);
    const colorScale = useRef();
    const composition = useRef();
    const selectedComposition = useRef();
    const title = useRef();
    const obsCat = searchTokens.filter(item => item.type === FEATURE_TYPE.OBS_CAT).map(item => item.id);
    const dimension = obsCat.length > 1 ? obsCat[obsCat.length - 1] : null;
    const obsCatString = JSON.stringify(obsCat);
    useEffect(() => {
        setNeedsUpdate(true);
    }, [categoricalNames, selection, obsCatString]);

    useEffect(() => {
        colorScale.current = getColorScale(embeddingData, dimension);
    }, [embeddingData, dimension]);


    useEffect(() => {
        if (needsUpdate && visible) {
            composition.current = getComposition(dataset, obsCat, cachedData, categoricalNames);
            selectedComposition.current = selection != null && selection.size > 0 ? getComposition(dataset, obsCat, cachedData, categoricalNames, selection) : null;
            title.current = dimension + ' composition in ' + obsCat.slice(0, obsCat.length - 1).join(', ');
            setNeedsUpdate(false);
        }
    }, [needsUpdate, visible, cachedData, categoricalNames, dataset, dimension, obsCat, selection]);

    if (dimension == null) {
        return null;
    }
    const textColor = darkMode ? 'white' : 'black';
    return <>{composition.current && <CompositionPlot seriesToValueToCounts={composition.current.seriesToValueToCounts}
                                                      dimension={dimension.current}
                                                      title={title.current}
                                                      colorScale={colorScale.current}
                                                      series={composition.current.series}
                                                      uniqueValues={composition.current.uniqueValues}
                                                      textColor={textColor}/>}

        {selectedComposition.current &&
        <CompositionPlot seriesToValueToCounts={selectedComposition.current.seriesToValueToCounts}
                         dimension={dimension}
                         title={title.current}
                         subtitle="selection"
                         colorScale={colorScale.current}
                         series={selectedComposition.current.series}
                         uniqueValues={selectedComposition.current.uniqueValues}
                         textColor={textColor}/>}
    </>;


}

const mapStateToProps = state => {
    return {
        cachedData: state.cachedData,
        categoricalNames: state.categoricalNames,
        darkMode: state.chartOptions.darkMode,
        dataset: state.dataset,
        embeddingData: state.embeddingData,
        searchTokens: state.searchTokens,
        selection: state.selection,
        visible: state.tab === 'composition'
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(CompositionPlots));


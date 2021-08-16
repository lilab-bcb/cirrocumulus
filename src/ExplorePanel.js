import AutocompleteVirtualized from "./AutocompleteVirtualized";
import {
    copyToClipboard,
    FEATURE_TYPE,
    getFeatureSets,
    NATSORT,
    SERVER_CAPABILITY_SAVE_FEATURE_SETS,
    splitSearchTokens
} from "./util";
import NumberIcon from "./NumberIcon";
import {InputLabel, Switch, Typography} from '@material-ui/core';
import FormControl from '@material-ui/core/FormControl';
import IconButton from '@material-ui/core/IconButton';
import Link from '@material-ui/core/Link';
import Tooltip from '@material-ui/core/Tooltip';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import FontDownloadRoundedIcon from '@material-ui/icons/FontDownloadRounded';
import {findIndex} from 'lodash';
import memoize from "memoize-one";
import React, {useState} from 'react';
import withStyles from '@material-ui/core/styles/withStyles';
import {
    deleteFeatureSet,
    deleteView,
    downloadSelectedIds,
    getDatasetFilterNames,
    getEmbeddingKey,
    getTraceKey,
    removeDatasetFilter,
    SAVE_FEATURE_SET_DIALOG,
    setActiveFeature,
    setCombineDatasetFilters,
    setDialog,
    setSearchTokens,
    setSelectedEmbedding,
    setTab,
    toggleEmbeddingLabel
} from './actions';
import {intFormat} from "./formatters";
import Chip from '@material-ui/core/Chip';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Divider from '@material-ui/core/Divider';
import Grid from '@material-ui/core/Grid';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import TextField from '@material-ui/core/TextField';
import CloudDownloadIcon from '@material-ui/icons/CloudDownload';
import HighlightOffIcon from '@material-ui/icons/HighlightOff';
import {connect} from 'react-redux';

const styles = theme => ({
    toolbar: {
        '& hr': {
            margin: theme.spacing(0, 0.5),
        }
    },
});
const getAnnotationOptions = memoize(
    (obs, obsCat) => {
        const options = [];
        obs.forEach(item => {
            options.push({group: 'Continuous', text: item, id: item});
        });
        obsCat.forEach(item => {
            options.push({group: 'Categorical', text: item, id: item});
        });
        options.sort((item1, item2) => {
            // const c = NATSORT(item1.group, item2.group);
            // if (c !== 0) {
            //     return c;
            // }
            return NATSORT(item1.text, item2.text);
        });
        return options;
    }
);
const getEmbeddingOptions = memoize(
    (embeddings) => {
        const options = [];
        embeddings.forEach(embedding => {
            options.push({
                text: embedding.name + (embedding.mode ? ' ' + embedding.mode : ''),
                id: getEmbeddingKey(embedding)
            });
        });

        options.sort((item1, item2) => {
            return NATSORT(item1.text.toLowerCase(), item2.text.toLowerCase());
        });
        return options;
    }
);
const getMetafeaturesOptions = memoize((items) => {
        if (items) {
            const options = items.slice();
            options.sort(NATSORT);
            return options;
        }
    }
);
const getFeatureSetOptions = memoize((items, categoricalNames) => {
        const options = items.map(item => ({group: item.category, text: item.name, id: item.id}));

        options.forEach(item => {
            let group = item.group;
            let index = group.lastIndexOf(' (');
            if (index !== -1) {
                group = group.substring(0, index);
            }
            let map = categoricalNames[group] || {};
            let newName = map[item.text];
            if (newName != null) {
                item.text = newName;
            }
        });
        options.sort((item1, item2) => {
            const c = NATSORT(item1.group, item2.group);
            if (c !== 0) {
                return c;
            }
            return NATSORT(item1.text, item2.text);
        });
        return options;
    }
);


function ExplorePanel(props) {


    const [featureSetAnchorEl, setFeatureSetAnchorEl] = useState(null);
    const [featureSet, setFeatureSet] = useState(null);
    const [featureSetView, setFeatureSetView] = useState(null);

    function onDatasetFilterChipDeleted(name) {
        props.removeDatasetFilter(name);
    }

    function onDatasetFilterCleared() {
        props.removeDatasetFilter(null);
    }

    function handleCombineDatasetFilters(event) {
        props.handleCombineDatasetFilters(event.target.checked ? 'or' : 'and');
    }

    function onFeaturesChange(event, value) {
        props.handleSearchTokens(value, FEATURE_TYPE.X);
    }

    function onMetafeaturesChange(event, value) {
        props.handleSearchTokens(value, FEATURE_TYPE.METAFEATURE);
    }

    function onObservationsIconClick(event, option) {
        props.handleEmbeddingLabel(option);
        event.stopPropagation();
    }

    function onObservationsChange(event, value) {
        let values = [];
        value.forEach(val => {
            if (val.text !== undefined) {
                values.push(val.text);
            } else {
                values.push(val);
            }
        });
        props.handleSearchTokens(values, FEATURE_TYPE.OBS);
    }

    function onSaveFeatureList() {
        props.handleDialog(SAVE_FEATURE_SET_DIALOG);
    }


    function onEmbeddingsChange(event, value) {
        const selection = [];
        const embeddingKeys = props.dataset.embeddings.map(item => getEmbeddingKey(item));
        value.forEach(val => {
            const id = val.id !== undefined ? val.id : val;
            const index = embeddingKeys.indexOf(id);
            let embedding = props.dataset.embeddings[index];
            selection.push(embedding);
        });
        props.handleEmbeddings(selection);
    }

    function onFeatureSetsChange(event, value) {
        let values = [];
        value.forEach(val => {
            if (val.id !== undefined) {
                values.push(val.id);
            } else {
                values.push(val);
            }
        });
        props.handleSearchTokens(values, FEATURE_TYPE.FEATURE_SET);
    }

    function onFeatureSetClick(event, option) {
        event.stopPropagation();
        const id = option.id;
        const target = event.target.closest(".MuiFormControl-root");
        let markers = props.markers;
        let newFeatureSet = null;
        for (let i = 0; i < markers.length; i++) {
            if (markers[i].id === id) {
                newFeatureSet = markers[i];
                break;
            }
        }
        if (newFeatureSet == null) {
            console.log(id + ' not found');
        }
        setFeatureSetAnchorEl(target);
        setFeatureSet(newFeatureSet);
    }

    function onFeatureSetMenuClose(event) {
        setFeatureSetAnchorEl(null);
        setFeatureSet(null);
    }

    function onViewFeatureSet(event) {
        event.stopPropagation();
        setFeatureSetAnchorEl(null);
        setFeatureSetView(true);
    }

    function onCloseViewFeatureSetDialog(event) {
        event.stopPropagation();
        setFeatureSet(null)
        setFeatureSetView(null)
    }

    function onDeleteFeatureSet(event) {
        event.stopPropagation();
        let searchTokens = props.searchTokens;
        const featureSetId = featureSet.id;
        let value = searchTokens.filter(token => token.type === FEATURE_TYPE.FEATURE_SET && token.value.id !== featureSetId);
        props.handleSearchTokens(value, FEATURE_TYPE.FEATURE_SET);
        props.handleDeleteFeatureSet(featureSetId);
        setFeatureSet(null)
        setFeatureSetAnchorEl(null);
    }

    function onFilterChipClicked(event) {
        onFeatureClick(event, event.target.innerText);
    }

    function onFeatureCopy(event) {
        event.preventDefault();
        event.stopPropagation();
        let searchTokens = props.searchTokens;
        copyToClipboard(searchTokens.filter(token => token.type === FEATURE_TYPE.X).map(item => item.value).join('\n'));
    }

    function onDownloadSelectedIds(event) {
        event.preventDefault();
        props.handleDownloadSelectedIds();
    }

    function onFeatureClick(event, option) {
        event.stopPropagation();
        const value = option.text !== undefined ? option.text : option;
        let galleryTraces = props.embeddingData.filter(traceInfo => traceInfo.active);
        for (let i = 0; i < galleryTraces.length; i++) {
            if (galleryTraces[i].name === value) {
                if (props.tab !== 'embedding') {
                    props.handleTab('embedding');
                }
                props.handleActiveFeature({
                    name: galleryTraces[i].name,
                    type: galleryTraces[i].featureType,
                    embeddingKey: getTraceKey(galleryTraces[i])
                });
                break;
            }
        }
    }

    const {
        categoricalNames,
        classes,
        combineDatasetFilters,
        dataset,
        datasetFilter,
        embeddingLabels,
        embeddings,
        markers,
        searchTokens,
        selection,
        serverInfo,
        tab,
    } = props;

    const datasetFilterKeys = getDatasetFilterNames(datasetFilter);
    datasetFilterKeys.sort(NATSORT);
    const splitTokens = splitSearchTokens(searchTokens);
    const featureSets = getFeatureSets(markers, splitTokens.featureSets);
    const featureOptions = dataset.features;
    const metafeatureOptions = getMetafeaturesOptions(dataset.metafeatures);
    const obsCat = dataset.obsCat;
    const obs = dataset.obs;
    const annotationOptions = getAnnotationOptions(obs, obsCat);
    const featureSetOptions = getFeatureSetOptions(markers, categoricalNames);
    const embeddingOptions = getEmbeddingOptions(dataset.embeddings);
    const selectedEmbeddings = getEmbeddingOptions(embeddings);
    return <>

        <Dialog
            open={Boolean(featureSetView)}
            onClose={onCloseViewFeatureSetDialog}
            aria-labelledby="view-set-dialog-title"
            fullWidth={true}
            maxWidth={'lg'}
        >
            <DialogTitle id="view-set-dialog-title">{featureSet ? featureSet.name : ''}</DialogTitle>
            <DialogContent>
                <TextField
                    value={featureSet ? featureSet.features.join('\n') : ''}
                    margin="dense"
                    fullWidth
                    readOnly={true}
                    variant="outlined"
                    multiline={true}
                />
            </DialogContent>
        </Dialog>
        <Menu
            id="feature-set-menu"
            anchorEl={featureSetAnchorEl}
            open={Boolean(featureSetAnchorEl)}
            onClose={onFeatureSetMenuClose}
        >
            <MenuItem onClick={onViewFeatureSet}>View</MenuItem>
            <MenuItem divider={true}/>
            <MenuItem disabled={featureSet && featureSet.readonly}
                      onClick={onDeleteFeatureSet}>Delete</MenuItem>

        </Menu>
        <div style={tab === 'embedding' || tab === 'distribution' || tab === 'composition' ? null : {display: 'none'}}>
            <div className={classes.section}>
                <Divider/>
                <Typography gutterBottom={false} component={"h1"}
                            style={{textTransform: 'uppercase', letterSpacing: '0.1em'}}>Explore</Typography>
                {tab === 'embedding' && embeddingOptions.length > 0 &&
                <FormControl className={classes.formControl}>
                    <AutocompleteVirtualized label={"Embeddings"}
                                             testId={'embeddings-input'}
                                             options={embeddingOptions}
                                             getChipTitle={(option) => {
                                                 return option.text;
                                             }}
                                             value={selectedEmbeddings}
                                             getChipText={(option) => option.text}
                                             getOptionLabel={(option) => option.text}
                                             getOptionSelected={(option, value) => findIndex(selectedEmbeddings, item => item.id === option.id) !== -1}
                                             onChange={onEmbeddingsChange}
                    />
                </FormControl>}
                {featureOptions.length > 0 && <FormControl className={classes.formControl}>
                    <AutocompleteVirtualized onChipClick={onFeatureClick}
                                             label={"Genes/Features"}
                                             testId={'genes-input'}
                                             options={featureOptions}
                                             value={splitTokens.X}
                                             onChange={onFeaturesChange}
                                             helperText={"Enter or paste list"}

                    />
                    <div><Link
                        style={{
                            float: 'right',
                            fontSize: '0.75rem',
                            transform: 'translateY(-50px)',
                            display: splitTokens.X.length === 0 ? 'none' : ''
                        }}
                        onClick={onFeatureCopy}>Copy</Link></div>

                </FormControl>}
                {annotationOptions.length > 0 && <FormControl className={classes.formControl}>
                    <AutocompleteVirtualized label={"Cell Metadata"}
                                             testId={'cell-meta-input'}
                                             options={annotationOptions}
                                             value={searchTokens.filter(token => token.type === FEATURE_TYPE.OBS_CAT || token.type === FEATURE_TYPE.OBS).map(token => token.value)}
                                             onChipClick={onFeatureClick}
                                             groupBy={false}
                                             getOptionLabel={(option) => option.text}
                                             getOptionIcon={(option) => option.group === 'Categorical' ?
                                                 <FontDownloadRoundedIcon style={{
                                                     marginRight: 2,
                                                     fontSize: '0.9rem'

                                                 }}/> : <NumberIcon style={{
                                                     marginRight: 2,
                                                     fontSize: '0.9rem'
                                                 }}/>}
                                             getChipIcon={(option) => {
                                                 return splitTokens.obsCat.indexOf(option) !== -1 ?
                                                     <FontDownloadRoundedIcon
                                                         onClick={(event) => {
                                                             onObservationsIconClick(event, option);
                                                         }}
                                                         title={"Toggle Show/Hide Labels"}
                                                         style={{
                                                             marginLeft: 4,
                                                             marginTop: 0,
                                                             marginRight: 0,
                                                             marginBottom: 0
                                                         }}
                                                         className={"MuiChip-deleteIcon MuiChip-deleteIconSmall" + (embeddingLabels.indexOf(option) !== -1 ? ' cirro-active' : '')}/> : null;
                                             }}
                                             getOptionSelected={(option, value) => option.id === value}
                                             onChange={onObservationsChange}/>
                </FormControl>}
                {metafeatureOptions && metafeatureOptions.length > 0 &&
                <FormControl className={classes.formControl}>
                    <AutocompleteVirtualized
                        label={"Metagenes/Features"}
                        options={metafeatureOptions}
                        value={splitTokens.metafeatures}
                        onChange={onMetafeaturesChange}
                    />
                </FormControl>}
                {<FormControl className={classes.formControl}>
                    <AutocompleteVirtualized label={"Sets"}
                                             testId={'sets-input'}
                                             options={featureSetOptions}
                                             value={featureSets}
                                             onChipClick={onFeatureSetClick}
                                             getChipTitle={(option) => {
                                                 return option.category + ', ' + option.name;
                                             }}
                                             getChipIcon={(option) => {
                                                 return <ArrowDropDownIcon onClick={(event) => {
                                                     onFeatureSetClick(event, option);
                                                 }}/>;
                                             }}
                                             groupBy={true}
                                             onChange={onFeatureSetsChange}
                                             getOptionSelected={(option, value) => option.id === value.id}
                                             getChipText={option => option.name}
                    />
                    {serverInfo.capabilities.has(SERVER_CAPABILITY_SAVE_FEATURE_SETS) && <div>
                        <Tooltip title={"Save Current Genes/Features"}>
                            <Link
                                style={{
                                    float: 'right',
                                    fontSize: '0.75rem',
                                    transform: 'translateY(-50px)',
                                    display: splitTokens.X.length === 0 ? 'none' : ''
                                }}
                                onClick={onSaveFeatureList}>Save</Link></Tooltip></div>}

                </FormControl>}
            </div>
            <div className={classes.section} style={{maxHeight: 500}}>
                <Divider inset='true'/>
                <Typography gutterBottom={false} component={"h1"}
                            style={{textTransform: 'uppercase'}}>Filters</Typography>
                <Grid component="label" alignContent={"flex-start"} container alignItems="center"
                      spacing={0}>
                    <Grid item><InputLabel shrink={true} variant={"standard"}>Combine</InputLabel></Grid>
                    <Grid item>AND</Grid>
                    <Grid item>
                        <Switch
                            size="small"
                            checked={combineDatasetFilters === 'or'}
                            onChange={handleCombineDatasetFilters}
                        />
                    </Grid>
                    <Grid item>OR</Grid>
                </Grid>
                {datasetFilterKeys.length > 0 && selection.size > 0 &&
                <>
                    <div style={{marginBottom: 2}}>
                        {intFormat(selection.size) + " / " + intFormat(dataset.shape[0]) + ": "}
                        {datasetFilterKeys.map(key => {
                            return <Chip
                                onDelete={() => {
                                    onDatasetFilterChipDeleted(key);
                                }}
                                onClick={onFilterChipClicked} size={"small"} variant={"default"}
                                style={{marginRight: 2, verticalAlign: 'bottom'}}
                                key={key}
                                label={key}

                            />;
                        })}
                        <Divider/>
                        <Grid container alignItems="center" className={classes.toolbar}
                              disabled={datasetFilterKeys.length === 0}>
                            <Tooltip title={"Clear All"}>
                                <IconButton size={'small'}
                                            onClick={onDatasetFilterCleared}><HighlightOffIcon/></IconButton>
                            </Tooltip>
                            <Tooltip title={"Download Selected IDs"}>
                                <IconButton size={'small'}
                                            onClick={onDownloadSelectedIds}><CloudDownloadIcon/></IconButton>
                            </Tooltip>
                        </Grid>
                    </div>
                </>
                }
            </div>
        </div>
    </>
}


const mapStateToProps = state => {
        return {
            activeFeature: state.activeFeature,
            categoricalNames: state.categoricalNames,
            combineDatasetFilters: state.combineDatasetFilters,
            dataset: state.dataset,
            datasetFilter: state.datasetFilter,
            datasetFilters: state.datasetFilters,
            embeddingData: state.embeddingData,
            embeddingLabels: state.embeddingLabels,
            embeddings: state.embeddings,
            markers: state.markers,
            savedDatasetFilter: state.savedDatasetFilter,
            searchTokens: state.searchTokens,
            selection: state.selection,
            serverInfo: state.serverInfo,
            tab: state.tab
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
    return {
        handleDialog: (value) => {
            dispatch(setDialog(value));
        },
        handleTab: (value) => {
            dispatch(setTab(value));
        },
        handleActiveFeature: (value) => {
            dispatch(setActiveFeature(value));
        },
        handleCombineDatasetFilters: (value) => {
            dispatch(setCombineDatasetFilters(value));
        },
        handleDownloadSelectedIds: () => {
            dispatch(downloadSelectedIds());
        },
        removeDatasetFilter: (filter) => {
            dispatch(removeDatasetFilter(filter));
        },
        handleEmbeddings: value => {
            dispatch(setSelectedEmbedding(value));
        },
        handleEmbeddingLabel: value => {
            dispatch(toggleEmbeddingLabel(value));
        },
        handleSearchTokens: (value, type) => {
            dispatch(setSearchTokens(value == null ? [] : value, type));
        },
        handleDeleteView: value => {
            dispatch(deleteView(value));
        },
        handleDeleteFeatureSet: value => {
            dispatch(deleteFeatureSet(value));
        }
    };
}

export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(ExplorePanel));
import {Tooltip} from '@material-ui/core';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import IconButton from '@material-ui/core/IconButton';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import TextField from '@material-ui/core/TextField';
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown';
import React from 'react';
import {intFormat, numberFormat} from './formatters';
import {NATSORT} from './util';

class CategoricalLegend extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {
            contextmenuEl: null,
            anchorEl: null,
            color: null,
            name: '',
            categoryValue: null,
            forceUpdate: false,
            menu: null
        };
    }

    handlePopoverClose = (e) => {
        this.setState({contextmenuEl: null, anchorEl: null});
    };

    handleColorChange = (e) => {
        this.setState({color: e.target.value});
    };
    handleNameChange = (e) => {
        this.setState({name: e.target.value});
    };

    handleColorChangeApply = (e) => {
        this.props.handleColorChange({
            name: this.props.name,
            value: this.state.categoryValue,
            color: this.state.color
        });
        this.setState({forceUpdate: !this.state.forceUpdate});
    };

    handleNameChangeApply = (e) => {
        this.props.handleNameChange({
            name: this.props.name,
            oldValue: this.state.categoryValue,
            value: this.state.name
        });
        this.setState({name: '', contextmenuEl: null, anchorEl: null});
    };

    handleEditColor = (e) => {
        this.setState((prevState, props) => ({
            anchorEl: prevState.contextmenuEl,
            contextmenuEl: null,
            menu: 'color'
        }), () => {
            // this.inputElement.click();
        });
    };
    handleEditName = (e) => {
        this.setState((prevState, props) => ({
            anchorEl: prevState.contextmenuEl,
            contextmenuEl: null,
            menu: 'name',
            name: ''
        }), () => {
            // this.inputElement.click();
        });
    };

    handleContextmenuClose = (e) => {
        this.setState({contextmenuEl: null});
    };

    handleClick = (value, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            this.props.handleClick({
                name: this.props.name,
                value: value,
                shiftKey: e.shiftKey,
                metaKey: e.ctrlKey || e.metaKey
            });
        }
    };

    handleContextmenu = (value, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            e.stopPropagation();
            this.setState({contextmenuEl: e.target, categoryValue: value, color: this.props.scale(value)});
        }
    };


    render() {
        const {
            scale,
            dataset,
            datasetFilter,
            name,
            featureSummary,
            maxHeight,
            globalFeatureSummary,
            nObs,
            nObsSelected,
            categoricalNames
        } = this.props;
        let clickEnabled = this.props.clickEnabled;
        const categoricalFilter = datasetFilter[name];
        const selectionSummary = featureSummary[name];
        let selectedDimensionToCount = {};
        if (selectionSummary != null) {
            for (let i = 0; i < selectionSummary.counts.length; i++) {
                selectedDimensionToCount[selectionSummary.categories[i]] = selectionSummary.counts[i];
            }
        }
        const globalDimensionSummary = globalFeatureSummary[name];
        const categories = globalDimensionSummary.categories.slice(0); // make a copy so that when sorting, counts stays in same order as categories
        const renamedCategories = categoricalNames[name] || {};
        if (dataset.categoryOrder && dataset.categoryOrder[name]) {
            const orderedCategories = dataset.categoryOrder[name];
            const categoryToIndex = new Map();
            for (let i = 0; i < orderedCategories.length; i++) {
                categoryToIndex.set(orderedCategories[i], i);
            }
            categories.sort((a, b) => categoryToIndex.get(a) - categoryToIndex.get(b));
        } else {
            categories.sort((a, b) => {
                let renamed1 = renamedCategories[a];
                if (renamed1 != null) {
                    a = renamed1;
                }
                let renamed2 = renamedCategories[b];
                if (renamed2 != null) {
                    b = renamed2;
                }
                return NATSORT(a, b);
            });
        }
        clickEnabled = clickEnabled && categories.length > 1;
        let style = {maxHeight: maxHeight, display: 'inline-block'};
        if (this.props.style) {
            style = Object.assign({}, style, this.props.style);
        }

        const renamedCategoryValue = renamedCategories[this.state.categoryValue] || this.state.categoryValue;
        return (
            <div className="cirro-chart-legend" style={style}>
                <Dialog open={Boolean(this.state.anchorEl)} onClose={this.handlePopoverClose}
                        aria-labelledby="edit-category-dialog-title">

                    {Boolean(this.state.anchorEl) && this.state.menu == 'color' && <>
                        <DialogTitle id="edit-category-dialog-title">Edit {renamedCategoryValue} Color</DialogTitle>
                        <DialogContent>
                            <input type="color" value={this.state.color}
                                   onChange={this.handleColorChange} style={{width: 100}}/>
                        </DialogContent>
                        <DialogActions>
                            <Button onClick={this.handlePopoverClose} color="primary">
                                Close
                            </Button>
                            <Button onClick={this.handleColorChangeApply} color="primary">
                                Apply
                            </Button>
                        </DialogActions>
                    </>}
                    {Boolean(this.state.anchorEl) && this.state.menu == 'name' && <>
                        <DialogTitle id="edit-category-dialog-title">Edit {renamedCategoryValue} Name</DialogTitle>
                        <DialogContent>
                            <TextField type="text" onChange={this.handleNameChange} value={this.state.name}/>
                        </DialogContent>
                        <DialogActions>
                            <Button onClick={this.handlePopoverClose} color="primary">
                                Cancel
                            </Button>
                            <Button onClick={this.handleNameChangeApply} color="primary">
                                OK
                            </Button>
                        </DialogActions>
                    </>}
                </Dialog>
                <Menu
                    anchorEl={this.state.contextmenuEl}
                    open={Boolean(this.state.contextmenuEl)}
                    onClose={this.handleContextmenuClose}
                >
                    <MenuItem onClick={this.handleEditName}>Edit Name</MenuItem>
                    <MenuItem onClick={this.handleEditColor}>Edit Color</MenuItem>
                </Menu>
                <table style={{textAlign: 'left', userSelect: 'none'}}>
                    <thead>
                    <tr>
                        {clickEnabled && <td></td>}
                        <td></td>
                        <td><small>{'all'}</small></td>
                        <td><small>{selectionSummary != null ? 'selection' : null}</small></td>
                    </tr>
                    </thead>
                    <tbody>
                    {categories.map((category) => {
                        const categoryIndex = globalDimensionSummary.categories.indexOf(category);
                        if (categoryIndex === -1) {
                            throw new Error(category + ' not found');
                        }
                        const opacity = categoricalFilter == null || categoricalFilter.value.indexOf(category) !== -1 ? 1 : 0.4;
                        const categoryCount = globalDimensionSummary.counts[categoryIndex];

                        const selectedCategoryCount = selectedDimensionToCount[category] || 0;

                        //            not-selected, selected
                        // in category      a       b
                        // not in category  c       d
                        // const a = categoryCount - selectedCategoryCount;
                        // const b = selectedCategoryCount;
                        // const c = nObs - nObsSelected - selectedCategoryCount;
                        // const d = nObsSelected - selectedCategoryCount;

                        const fractionSelected = selectionSummary == null ? 0 : selectedCategoryCount / nObsSelected;
                        const globalTitle = numberFormat(100 * categoryCount / nObs) + '%';
                        let categoryText = category;
                        let renamed = renamedCategories[category];
                        if (renamed !== undefined) {
                            categoryText = renamed;
                        }
                        const selectionTitle = selectionSummary == null ? null : numberFormat(100 * fractionSelected) + '%';
                        return <tr
                            style={{cursor: clickEnabled ? 'pointer' : null, opacity: opacity}}
                            onContextMenu={(e) => this.handleContextmenu(category, e)}
                            onClick={(e) => this.handleClick(category, e)} key={category}>
                            {clickEnabled && <td>
                                <div style={{
                                    display: 'inline-block',
                                    width: 10,
                                    height: 10,
                                    background: scale(category)
                                }}/>

                            </td>}
                            <td>
                                <div style={{
                                    maxWidth: 140,
                                    textOverflow: 'ellipsis',
                                    overflow: 'hidden',
                                    display: 'inline-block',
                                    userSelect: 'none'
                                }} title={'' + categoryText}>{'' + categoryText}</div>
                                <IconButton style={{padding: 0, fontSize: 14}} size="small"
                                            onClick={(e) => this.handleContextmenu(category, e)} aria-label="Menu"
                                            aria-haspopup="true">
                                    <ArrowDropDownIcon fontSize={"inherit"}/>
                                </IconButton>
                            </td>

                            <td>
                                <Tooltip title={globalTitle}>
                                    <span>{intFormat(categoryCount)}</span>
                                </Tooltip>

                            </td>
                            {selectionSummary && <td>
                                <Tooltip title={selectionTitle}>
                                    <div>{intFormat(selectedCategoryCount)}</div>
                                </Tooltip>
                            </td>}
                        </tr>;
                    })
                    }</tbody>

                </table>
            </div>
        );
    }
}

export default CategoricalLegend;



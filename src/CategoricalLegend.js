import {Tooltip} from '@material-ui/core';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';

import {scaleLinear} from 'd3-scale';
import React from 'react';
import {intFormat, numberFormat} from './formatters';

class CategoricalLegend extends React.PureComponent {

    constructor(props) {
        super(props);
        this.state = {contextmenuEl: null, anchorEl: null, color: null, colorValue: null, forceUpdate: false};
    }


    handlePopoverClose = (e) => {
        this.setState({contextmenuEl: null, anchorEl: null});
    };


    handleColorChange = (e) => {
        this.setState({color: e.target.value});
    };

    handleColorChangeApply = (e) => {
        this.props.handleColorChange({name: this.props.name, value: this.state.colorValue, color: this.state.color});
        this.setState({forceUpdate: !this.state.forceUpdate});
    };

    handleEditColor = (e) => {
        this.setState((prevState, props) => ({
            anchorEl: prevState.contextmenuEl,
            contextmenuEl: null
        }), () => {
            // this.inputElement.click();
        });
    };

    handleContextmenuClose = (e) => {
        this.setState({contextmenuEl: null});
    };

    handleClick = (value, index, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            this.props.handleClick({name: this.props.name, value: value, shiftKey: e.shiftKey, metaKey: e.metaKey});
        }
    };

    handleContextmenu = (value, index, e) => {
        if (this.props.clickEnabled) {
            e.preventDefault();
            this.setState({contextmenuEl: e.target, colorValue: value, color: this.props.scale(value)});
        }
    };

    draw(context) {
        const {scale, name, globalFeatureSummary} = this.props;
        const globalDimensionSummary = globalFeatureSummary[name];
        const categories = globalDimensionSummary.categories;
        context.font = '14px Helvetica';
        context.textAlign = 'left';
        for (let i = 0; i < categories.length; i++) {
            const category = categories[i];
            context.fillStyle = scale(category);
            context.fillRect(0, i * 12, 10, 10);
            context.fillStyle = 'black';
            context.fillText('' + category, 12, i * 12);
        }
    }

    getSize(context) {
        const {name, globalFeatureSummary} = this.props;
        context.font = '14px Helvetica';
        const globalDimensionSummary = globalFeatureSummary[name];
        const categories = globalDimensionSummary.categories;
        let maxWidth = context.measureText(name).width;
        categories.forEach(value => maxWidth = Math.max(maxWidth, context.measureText(value).width));
        return {width: maxWidth + 14, height: categories.length * 12};
    }

    render() {
        const {scale, datasetFilter, name, featureSummary, maxHeight, globalFeatureSummary, nObs, nObsSelected} = this.props;
        let clickEnabled = this.props.clickEnabled;
        const categoricalFilterValues = datasetFilter[name];
        const selectionSummary = featureSummary[name];
        let selectedDimensionToCount = {};
        if (selectionSummary != null) {
            for (let i = 0; i < selectionSummary.counts.length; i++) {
                selectedDimensionToCount[selectionSummary.categories[i]] = selectionSummary.counts[i];
            }
        }
        const globalDimensionSummary = globalFeatureSummary[name];
        // if (globalDimensionSummary.max == null) {
        //     let max = 0;
        //     let min = Number.MAX_VALUE;
        //     for (let i = 0; i < globalDimensionSummary.counts.length; i++) {
        //         max = Math.max(max, globalDimensionSummary[i]);
        //         min = Math.min(min, globalDimensionSummary[i]);
        //     }
        //     globalDimensionSummary.min = min;
        //     globalDimensionSummary.max = max;
        // }
        const categories = globalDimensionSummary.categories;
        clickEnabled = clickEnabled && categories.length > 1;

        let maxSize = 60;
        const fractionScale = scaleLinear().domain([0, 1]).range([0, maxSize]).clamp(true);
        return (
            <div className="cirro-chart-legend" style={{
                maxHeight: maxHeight
            }}>
                <Dialog open={Boolean(this.state.anchorEl)} onClose={this.handlePopoverClose}
                        aria-labelledby="edit-category-color-dialog-title">
                    <DialogTitle id="edit-category-color-dialog-title">Edit {this.state.colorValue} Color</DialogTitle>
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
                </Dialog>
                <b>{name}</b> <small>({categories.length})</small>


                {/*<Popover*/}
                {/*    open={Boolean(this.state.anchorEl)}*/}
                {/*    anchorEl={this.state.anchorEl}*/}
                {/*    onClose={this.handlePopoverClose}*/}
                {/*    anchorOrigin={{*/}
                {/*        vertical: 'bottom',*/}
                {/*        horizontal: 'center',*/}
                {/*    }}*/}
                {/*    transformOrigin={{*/}
                {/*        vertical: 'top',*/}
                {/*        horizontal: 'center',*/}
                {/*    }}*/}
                {/*>*/}
                {/*    <Typography>Edit {this.state.colorValue} Color</Typography>*/}
                {/*    <input ref={input => this.inputElement = input} type="color" value={this.state.color}*/}
                {/*           onChange={this.handleColorChange} style={{width: 100}}/>*/}
                {/*    <Button onClick={handleClose} color="primary">*/}
                {/*        Cancel*/}
                {/*    </Button>*/}
                {/*    <Button onClick={handleClose} color="primary">*/}
                {/*        Apply*/}
                {/*    </Button>*/}
                {/*</Popover>*/}
                <Menu
                    anchorEl={this.state.contextmenuEl}
                    open={Boolean(this.state.contextmenuEl)}
                    onClose={this.handleContextmenuClose}
                >
                    <MenuItem onClick={this.handleEditColor}>Edit Color</MenuItem>
                </Menu>
                <table>
                    <thead>
                    <tr>
                        {clickEnabled && <td></td>}
                        <td></td>
                        <td><small>{'all'}</small></td>
                        <td><small>{selectionSummary != null ? 'selection' : null}</small></td>
                    </tr>
                    </thead>
                    <tbody>
                    {categories.map((category, i) => {

                        const opacity = categoricalFilterValues == null || categoricalFilterValues.indexOf(category) !== -1 ? 1 : 0.4;
                        // const fractionUnselected = unselectedCounts != null ? unselectedCounts[i] / unselectedTotal : null;
                        // const unselectedSize = unselectedCounts == null ? 0 : fractionScale(fractionUnselected);
                        // const unselectedTitle = unselectedCounts == null ? null : intFormat(unselectedCounts[i]) + ' / ' + intFormat(unselectedTotal) + (unselectedCounts[i] > 0 ? (' ( ' + numberFormat(100 * fractionUnselected) + '%)') : '');
                        const count = globalDimensionSummary.counts[i];
                        const selectedCount = selectedDimensionToCount[category] || 0;
                        const fractionSelected = selectionSummary == null ? 0 : selectedCount / nObsSelected;
                        // const selectedSize = fractionScale(fractionSelected);
                        // const globalSize = fractionScale(count / nObs);
                        const globalTitle = numberFormat(100 * count / nObs) + '%';
                        const selectionTitle = selectionSummary == null ? null : numberFormat(100 * fractionSelected) + '%';
                        return <tr
                            style={{cursor: clickEnabled ? 'pointer' : null, opacity: opacity}}
                            onContextMenu={(e) => this.handleContextmenu(category, i, e)}
                            onClick={(e) => this.handleClick(category, i, e)} key={category}>
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
                                }} title={'' + category}>{'' + category}</div>
                            </td>

                            <td>
                                <Tooltip title={globalTitle}>
                                    <div>{intFormat(count)}</div>
                                </Tooltip>
                                {/*<div*/}
                                {/*    title={globalTitle}*/}
                                {/*    style={{*/}
                                {/*        display: 'inline-block',*/}
                                {/*        position: 'relative',*/}
                                {/*        width: maxSize,*/}
                                {/*        border: '1px solid black',*/}
                                {/*        height: 9*/}
                                {/*    }}>*/}

                                {/*    <div style={{*/}
                                {/*        position: 'absolute',*/}
                                {/*        width: globalSize,*/}
                                {/*        left: 0,*/}
                                {/*        top: 0,*/}
                                {/*        backgroundColor: 'LightGrey',*/}
                                {/*        height: 9*/}
                                {/*    }}/>*/}
                                {/*</div>*/}
                            </td>
                            {selectionSummary && <td>
                                <Tooltip title={selectionTitle}>
                                    <div>{intFormat(selectedCount)}</div>
                                </Tooltip>
                                {/*<div*/}
                                {/*    title={selectionTitle}*/}
                                {/*    style={{*/}
                                {/*        display: 'inline-block',*/}
                                {/*        position: 'relative',*/}
                                {/*        width: maxSize,*/}
                                {/*        border: '1px solid black',*/}
                                {/*        height: 9*/}
                                {/*    }}>*/}

                                {/*    <div style={{*/}
                                {/*        position: 'absolute',*/}
                                {/*        width: selectedSize,*/}
                                {/*        left: 0,*/}
                                {/*        top: 0,*/}
                                {/*        backgroundColor: 'LightGrey',*/}
                                {/*        height: 9*/}
                                {/*    }}/>*/}
                                {/*</div>*/}
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



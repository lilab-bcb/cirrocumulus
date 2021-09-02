import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import TextField from '@material-ui/core/TextField';
import React, {useState} from 'react';
import {NATSORT} from './util';
import {IconButton, ListItem, ListItemText} from '@material-ui/core';
import {intFormat} from './formatters';
import {FixedSizeList} from 'react-window';
import MenuIcon from '@material-ui/icons/Menu';


export default function CategoricalLegend(props) {
    const [contextMenu, setContextMenu] = useState(null);
    const [tmpName, setTmpName] = useState('');
    const [menu, setMenu] = useState(null);
    const [color, setColor] = useState(null);
    const [sort, setSort] = useState(null);
    const [originalCategory, setOriginalCategory] = useState(null);
    const [forceUpdate, setForceUpdate] = useState(false);

    function handleDialogClose(e) {
        setMenu(null);
    }

    function handleColorChange(e) {
        setColor(e.target.value);
    }

    function handleNameChange(e) {
        setTmpName(e.target.value);
    }

    function handleColorChangeApply(e) {
        props.handleColorChange({
            name: props.name,
            value: originalCategory,
            color: color
        });
        setForceUpdate(!forceUpdate);
    }

    function handleNameChangeApply(e) {
        props.handleNameChange({
            name: props.name,
            originalValue: originalCategory,
            newValue: tmpName
        });
        setMenu(null);
        setTmpName('');
    }

    function handleEditColor(e) {
        setContextMenu(null);
        setMenu('color');
    }

    function handleEditName(e) {
        setContextMenu(null);
        setTmpName('');
        setMenu('name');

    }

    function handleContextmenuClose(e) {
        setContextMenu(null);
    }

    function onRowClick(event, category) {
        event.preventDefault();
        props.handleClick({
            name: props.name,
            value: category,
            shiftKey: event.shiftKey,
            metaKey: event.ctrlKey || event.metaKey
        });
    }

    function onContextmenu(event, originalCategory) {
        event.preventDefault();
        event.stopPropagation();

        setContextMenu({
            mouseX: event.clientX - 2,
            mouseY: event.clientY - 4
        });
        setOriginalCategory(originalCategory);
        setColor(props.scale(originalCategory));
    }


    const {
        scale,
        datasetFilter,
        name,
        height,
        featureSummary,
        globalFeatureSummary,
        nObs,
        nObsSelected,
        categoricalNames
    } = props;

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
    if (sort === 'ascending' || sort === 'descending') {
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
        if (sort === 'descending') {
            categories.reverse();
        }
    }

    //
    // function cellRenderer({
    //                           cellData,
    //                           columnData,
    //                           columnIndex,
    //                           dataKey,
    //                           isScrolling,
    //                           rowData,
    //                           rowIndex
    //                       }) {
    //     return (
    //         <TableCell
    //             component="div"
    //             className={clsx(classes.tableCell, classes.flexContainer, {
    //                 [classes.noClick]: onRowClick == null
    //             })}
    //             variant="body"
    //             style={{height: rowHeight}}
    //             align={'left'}
    //         >
    //             {columnIndex === 0 && <div style={{
    //                 display: 'inline-block',
    //                 width: 11,
    //                 height: 11,
    //                 background: scale(rowData.category)
    //             }}></div>}
    //             {columnIndex === 1 && <div>
    //                 <div style={{
    //                     maxWidth: maxCategoryWidth,
    //                     textOverflow: 'ellipsis',
    //                     overflow: 'hidden',
    //                     display: 'inline-block',
    //                     userSelect: 'none'
    //                 }} title={'' + rowData.renamedCategory}>{'' + rowData.renamedCategory}</div>
    //                 <IconButton style={{padding: 0, fontSize: 14}} size="small"
    //                             onClick={(e) => handleContextmenu(rowData.category, rowData.renamedCategory, e)}
    //                             aria-label="Menu"
    //                             aria-haspopup="true">
    //                     <ArrowDropDownIcon fontSize={"inherit"}/>
    //                 </IconButton></div>}
    //
    //             {columnIndex === 2 &&
    //             <span>{(selectionSummary == null ? '' : intFormat(selectedDimensionToCount[rowData.category] || 0) + ' / ') + intFormat(globalDimensionSummary.counts[rowData.categoryIndex])}</span>}
    //         </TableCell>
    //     );
    // }


    function renderRow(props) {
        const {index, style} = props;
        const category = categories[index];
        const categoryIndex = globalDimensionSummary.categories.indexOf(category);
        const isSelected = categoricalFilter != null && categoricalFilter.value.indexOf(category) !== -1;
        let renamedCategory = renamedCategories[category];
        if (renamedCategory == null) {
            renamedCategory = category;
        }
        return (
            <ListItem disableGutters={true} divider={true} dense={true}
                      onContextMenu={event => onContextmenu(event, category)}
                      onClick={event => onRowClick(event, category)}
                      selected={isSelected} button style={style}
                      key={index}>

                <div style={{
                    marginBottom: 15,
                    marginRight: 2,
                    display: 'inline-block',
                    width: 12,
                    height: 12,
                    background: scale(category)
                }}></div>

                <ListItemText title={renamedCategory} primaryTypographyProps={{noWrap: true}} primary={renamedCategory}
                              secondary={(selectionSummary == null ? '' : intFormat(selectedDimensionToCount[category] || 0) + ' / ') + intFormat(globalDimensionSummary.counts[categoryIndex])}/>
                <IconButton onClick={event => onContextmenu(event, category)} edge="end" aria-label="menu">
                    <MenuIcon></MenuIcon>
                </IconButton>
            </ListItem>
        );
    }


    const renamedCategoryValue = renamedCategories[originalCategory] || originalCategory;
    if (!props.visible) {
        return null;
    }

    return (
        <>
            <div data-testid="categorical-legend">
                <FixedSizeList height={height} width={250} itemSize={40}
                               itemCount={categories.length}>
                    {renderRow}
                </FixedSizeList>
            </div>

            <Dialog open={Boolean(menu)} onClose={handleDialogClose}
                    aria-labelledby="edit-category-dialog-title">
                {menu == 'color' && <>
                    <DialogTitle id="edit-category-dialog-title">Edit {renamedCategoryValue} Color</DialogTitle>
                    <DialogContent>
                        <input type="color" value={color}
                               onChange={handleColorChange} style={{width: 100}}/>
                    </DialogContent>
                    <DialogActions>
                        <Button onClick={handleDialogClose} color="primary">
                            Close
                        </Button>
                        <Button onClick={handleColorChangeApply} color="primary">
                            Apply
                        </Button>
                    </DialogActions>
                </>}
                {menu == 'name' && <>
                    <DialogTitle id="edit-category-dialog-title">Edit {renamedCategoryValue} Name</DialogTitle>
                    <DialogContent>
                        <TextField type="text" onChange={handleNameChange} value={tmpName}/>
                    </DialogContent>
                    <DialogActions>
                        <Button onClick={handleDialogClose} color="primary">
                            Cancel
                        </Button>
                        <Button onClick={handleNameChangeApply} color="primary">
                            OK
                        </Button>
                    </DialogActions>
                </>}
            </Dialog>
            <Menu
                anchorReference="anchorPosition"
                anchorPosition={
                    contextMenu != null
                        ? {top: contextMenu.mouseY, left: contextMenu.mouseX}
                        : undefined
                }
                open={Boolean(contextMenu)}
                onClose={handleContextmenuClose}
            >
                <MenuItem onClick={handleEditName}>Edit Name</MenuItem>
                <MenuItem onClick={handleEditColor}>Edit Color</MenuItem>
            </Menu>
        </>
    );

}


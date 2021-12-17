import {Typography} from '@mui/material';
import Chip from '@mui/material/Chip';
import {styled, useTheme} from '@mui/material/styles';
import TextField from '@mui/material/TextField';
import useMediaQuery from '@mui/material/useMediaQuery';
import Autocomplete, {autocompleteClasses} from '@mui/material/Autocomplete';
import PropTypes from 'prop-types';
import React from 'react';
import {sortableContainer, sortableElement} from 'react-sortable-hoc';
import {VariableSizeList} from 'react-window';
import Popper from '@mui/material/Popper';
import Link from '@mui/material/Link';


const LISTBOX_PADDING = 0; // px

function getTextMatch(text, inputValue) {
    inputValue = inputValue.toLowerCase();
    if (inputValue !== '') {
        const index = text.toLowerCase().indexOf(inputValue);
        if (index !== -1) { // bold the matching text when searching
            const inputLength = inputValue.length;
            const before = text.substring(0, index);
            const match = text.substring(index, index + inputLength);
            const after = text.substring(index + inputLength);
            return [before, match, after];
        }
    }
}

function renderRow(props) {
    const {data, index, style} = props;
    const dataSet = data[index];
    const inlineStyle = {
        ...style,
        top: style.top + LISTBOX_PADDING
    };


    if (dataSet.hasOwnProperty('group')) {
        if (dataSet.group === '') {
            return null;
        }
        inlineStyle.whiteSpace = 'nowrap';
        const textMatch = dataSet.inputValue != null ? getTextMatch(dataSet.group, dataSet.inputValue) : null;
        if (textMatch) {
            return <Typography component="li" {...dataSet[0]} style={inlineStyle}>
                {textMatch[0]}<b>{textMatch[1]}</b>{textMatch[2]} {dataSet.selectGroup && <Link href="#"
                                                                                                onClick={e => dataSet.selectGroup(e, dataSet.group)}>All</Link>}
            </Typography>;
        }
        return <Typography component="li" {...dataSet[0]} style={inlineStyle}>
            {dataSet.group} {dataSet.selectGroup &&
            <Link href="#" onClick={e => dataSet.selectGroup(e, dataSet.group)}>All</Link>}
        </Typography>;

    }

    const item = dataSet[1];
    const text = item.text != null ? item.text : item;
    const icon = item.icon != null ? item.icon : null;
    const inputValue = dataSet[2];
    const textMatch = getTextMatch(text, inputValue);
    if (textMatch) {
        return <Typography component="li" {...dataSet[0]} title={text} noWrap style={inlineStyle}>
            {icon}{textMatch[0]}<b>{textMatch[1]}</b>{textMatch[2]}
        </Typography>;
    }
    return <Typography component="li" {...dataSet[0]} title={text} noWrap style={inlineStyle}>
        {icon}{text}
    </Typography>;

}

const OuterElementContext = React.createContext({});

const OuterElementType = React.forwardRef((props, ref) => {
    const outerProps = React.useContext(OuterElementContext);
    return <div ref={ref} {...props} {...outerProps} />;
});

function useResetCache(data) {
    const ref = React.useRef(null);
    React.useEffect(() => {
        if (ref.current != null) {
            ref.current.resetAfterIndex(0, true);
        }
    }, [data]);
    return ref;
}

// Adapter for react-window
const ListboxComponent = React.forwardRef(function ListboxComponent(props, ref) {
    const {children, ...other} = props;
    let itemData = [];

    children.forEach((item) => {
        itemData.push(item);
        if (item.children) {
            itemData = itemData.concat(item.children);
        }
    });

    const theme = useTheme();
    const smUp = useMediaQuery(theme.breakpoints.up('sm'), {
        noSsr: true
    });

    const itemCount = itemData.length;
    const itemSize = smUp ? 24 : 30;

    const getChildSize = (child) => {
        if (child.hasOwnProperty('group')) {
            return child.group === '' ? 0 : 30;
        }

        return itemSize;
    };

    const getHeight = () => {
        if (itemCount > 8) {
            return 8 * itemSize;
        }
        return itemData.map(getChildSize).reduce((a, b) => a + b, 0);
    };

    const gridRef = useResetCache(itemCount);

    return (
        <div ref={ref}>
            <OuterElementContext.Provider value={other}>
                <VariableSizeList
                    itemData={itemData}
                    height={getHeight() + 2 * LISTBOX_PADDING}
                    width="100%"
                    ref={gridRef}
                    outerElementType={OuterElementType}
                    innerElementType="ul"
                    itemSize={(index) => getChildSize(itemData[index])}
                    overscanCount={5}
                    itemCount={itemCount}
                >
                    {renderRow}
                </VariableSizeList>
            </OuterElementContext.Provider>
        </div>
    );
});

ListboxComponent.propTypes = {
    children: PropTypes.node
};

const StyledPopper = styled(Popper)({
    [`& .${autocompleteClasses.listbox}`]: {
        boxSizing: 'border-box',
        '& ul': {
            padding: 0,
            margin: 0
        }
    }
});


export default function AutocompleteVirtualized(props) {
    const ref = React.createRef();
    const inputValueRef = React.createRef('');
    const {value, options, onChange} = props;

    function onDrop(event) {
        event.preventDefault();
        event.stopPropagation();
        const files = event.dataTransfer.files;
        const reader = new FileReader();
        reader.onload = function (event) {
            let text = event.target.result.trim();
            let tokens = text.split(/[\n,\t]/);
            enterTokens(event, tokens);
        };

        reader.onerror = function (event) {
            alert("Unable to read file.");
        };

        reader.readAsText(files[0]);
        showDragIndicator(false);
    }

    function enterTokens(event, tokens) {
        const results = value;
        tokens = tokens.map(token => token.toLowerCase().trim().replace(/["']/g, '')).filter(token => token !== '');
        if (tokens.length > 0) {
            let textToOption = new Map();
            for (let i = 0; i < options.length; i++) {
                const option = options[i];
                const text = option.text != null ? option.text : option;
                const textLowerCase = text.toLowerCase();
                textToOption.set(textLowerCase, option);
            }
            results.forEach(option => { // remove existing
                const text = option.text != null ? option.text : option;
                textToOption.delete(text.toLowerCase());
            });
            tokens.forEach(token => {
                const option = textToOption.get(token);
                if (option != null) {
                    textToOption.delete(token); // delete so that we don't add new token 2x
                    results.push(option);
                }
            });
        }
        onChange(event, results);
    }

    function onPaste(event) {
        let text = event.clipboardData.getData('text/plain');
        if (text != null && text.length > 0) {
            event.preventDefault();
            event.stopPropagation();
            text = text.trim();
            let tokens = text.split(/[\n,\t]/);
            enterTokens(event, tokens);
        }
    }

    function showDragIndicator(show) {
        ref.current.style.background = show ? '#1976d2' : '';
    }

    function onDragOver(event) {
        if (event.dataTransfer.items.length > 0) {
            event.preventDefault();
            event.stopPropagation();
            showDragIndicator(true);
        }

    }

    function onDragEnd(event) {
        event.preventDefault();
        event.stopPropagation();
        showDragIndicator(false);
    }


    const onChipClick = props.onChipClick ? (event, option) => {
        props.onChipClick(event, option);
    } : null;

    let getOptionSelected = props.getOptionSelected;
    if (getOptionSelected == null) {
        getOptionSelected = props.groupBy ? (option, value) => option.id === value.id && option.group === value.group : (option, value) => option === value;
    }
    let getChipText = props.getChipText;
    if (getChipText == null) {
        getChipText = (option) => option;
    }
    let getChipIcon = props.getChipIcon;
    if (getChipIcon == null) {
        getChipIcon = (option) => null;
    }


    let getChipTitle = props.getChipTitle;
    if (getChipTitle == null) {
        getChipTitle = (option) => null;
    }
    let getOptionLabel = props.getOptionLabel;
    if (getOptionLabel == null) {
        getOptionLabel = props.groupBy ? (option) => option.text : (option) => option;
    }

    const selectGroup = props.selectGroup ? (event, group) => {
        const uniqueValues = new Set(value);
        options.forEach(option => {
            if (option.group === group && !uniqueValues.has(option.id)) {
                uniqueValues.add(option.id);
            }
        });
        // ensure everything in group is added
        // newValue.splice(index, 1);
        props.onChange(event, Array.from(uniqueValues));
    } : null;

    const filterOptions = (options, {inputValue}) => {
        inputValue = inputValue.trim().toLowerCase();
        if (inputValue === '') {
            return options;
        }
        const inputLength = inputValue.length;
        const exactMatches = [];
        const startsWithMatches = [];
        const containsMatches = [];
        for (let i = 0, n = options.length; i < n; i++) {
            const option = options[i];
            const text = getOptionLabel(option).toLowerCase();
            const index = text.indexOf(inputValue);
            if (index === 0) {
                if (text.length === inputLength) {
                    exactMatches.push(option);
                } else {
                    startsWithMatches.push(option);
                }
            } else if (index !== -1) {
                containsMatches.push(option);
            } else if (props.groupBy) { // check group
                const groupIndex = option.group.toLowerCase().indexOf(inputValue);
                if (groupIndex !== -1) {
                    exactMatches.push(option);
                }
            }
        }
        return exactMatches.concat(startsWithMatches).concat(containsMatches);
    };
    const handleTagDelete = (event, index) => {
        const newValue = value.slice();
        newValue.splice(index, 1);
        props.onChange(event, newValue);
    };

    const onSortEnd = (event) => {
        const newValue = value.slice();
        const tmp = newValue[event.newIndex];
        newValue[event.newIndex] = newValue[event.oldIndex];
        newValue[event.oldIndex] = tmp;
        onChange(event, newValue);
    };

    const SortableItem = sortableElement(({option, sortIndex}) => {
        return (
            <Chip
                tabIndex="-1"
                key={sortIndex}
                style={{zIndex: 1000000, maxWidth: 216, overflow: 'hidden', textOverflow: 'ellipsis'}}
                onDelete={event => handleTagDelete(event, sortIndex)}
                onClick={onChipClick ? event => onChipClick(event, option) : null}
                label={getChipText(option)}
                title={getChipTitle(option)}
                size="small"
                icon={getChipIcon(option)}/>
        );
    });

    const SortableList = sortableContainer(({items}) => {
        return (
            <ul style={{padding: 0, marginTop: 0, marginBottom: 0, maxWidth: 240, overflow: 'hidden'}}>
                {items.map((option, index) => {
                        return <SortableItem key={index} option={option} index={index} sortIndex={index}/>;
                    }
                )}
            </ul>
        );
    });
    const multiple = props.multiple != null ? props.multiple : true;

    return <>
        <Autocomplete
            data-testid={props.testId}
            multiple={multiple}
            ref={ref}
            size={"small"}
            disableListWrap
            freeSolo={props.freeSolo != null ? props.freeSolo : false}
            blurOnSelect={true}
            openOnFocus={false}
            autoHighlight={true}
            filterOptions={filterOptions}
            isOptionEqualToValue={getOptionSelected}
            value={value}
            filterSelectedOptions={true}
            getOptionLabel={getOptionLabel}
            groupBy={props.groupBy}
            ChipProps={{size: 'small'}}
            ListboxComponent={ListboxComponent}
            PopperComponent={StyledPopper}
            options={options}
            onChange={onChange}
            renderTags={(value, getTagProps) =>
                null
            }
            renderInput={(params) => (
                <TextField {...params} label={props.label} helperText={props.helperText} fullWidth={true}
                           margin={"dense"}/>
            )}
            renderGroup={(params) => {
                params.selectGroup = selectGroup;
                params.inputValue = inputValueRef.current;
                return params;
            }}
            onInputChange={(event, value, reason) => {
                inputValueRef.current = value;
            }}
            renderOption={(props, option, {inputValue}) => [props, option, inputValue]}
            onPaste={onPaste}
            onDrop={onDrop}
            onDragOver={onDragOver}
            onDragEnd={onDragEnd}
            onDragLeave={onDragEnd}
        />

        {multiple && <SortableList
            distance={2}
            onSortEnd={onSortEnd}
            axis="xy" items={props.value}
        />}

    </>;
};



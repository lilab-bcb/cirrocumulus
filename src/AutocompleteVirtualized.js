import {Typography} from '@material-ui/core';
import Chip from '@material-ui/core/Chip';
import ListSubheader from '@material-ui/core/ListSubheader';
import {makeStyles, useTheme} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import Autocomplete from '@material-ui/lab/Autocomplete';
import PropTypes from 'prop-types';
import React from 'react';
import {VariableSizeList} from 'react-window';


const LISTBOX_PADDING = 0; // px

function renderRow(props) {
    const {data, index, style} = props;
    return React.cloneElement(data[index], {
        style: {
            ...style,
            top: style.top + LISTBOX_PADDING,
        },
    });
}

const OuterElementContext = React.createContext({});

const OuterElementType = React.forwardRef((props, ref) => {
    const outerProps = React.useContext(OuterElementContext);
    return <div ref={ref} {...props} {...outerProps} />;
});

// Adapter for react-window
const ListboxComponent = React.forwardRef(function ListboxComponent(props, ref) {
    const {children, ...other} = props;
    const itemData = React.Children.toArray(children);
    const theme = useTheme();
    const smUp = useMediaQuery(theme.breakpoints.up('sm'), {noSsr: true});
    const itemCount = itemData.length;
    const itemSize = smUp ? 24 : 36;

    const getChildSize = (child) => {
        // if (React.isValidElement(child) && child.type === ListSubheader) {
        //     return 48;
        // }

        return itemSize;
    };

    const getHeight = () => {
        if (itemCount > 8) {
            return 8 * itemSize;
        }
        return itemData.map(getChildSize).reduce((a, b) => a + b, 0);
    };

    return (
        <div ref={ref}>
            <OuterElementContext.Provider value={other}>
                <VariableSizeList
                    itemData={itemData}
                    height={getHeight() + LISTBOX_PADDING}
                    width="100%"
                    key={itemCount}
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
    children: PropTypes.node,
};


const useStyles = makeStyles({
    paper: {width: 230},
    listbox: {
        '& ul': {
            padding: 0,
            margin: 0,
        },
    },
});


const renderGroup = (params) => [
    <ListSubheader disableGutters component="div">
        <Typography noWrap>{params.group}</Typography>
    </ListSubheader>,
    params.children,
];


export default function AutocompleteVirtualized(props) {
    const classes = useStyles();
    const ref = React.createRef();

    function onDrop(event) {
        event.preventDefault();
        event.stopPropagation();
        let dt = event.dataTransfer;
        let files = dt.files;
        let reader = new FileReader();
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
    };

    function enterTokens(event, tokens) {
        const results = props.value;
        tokens = tokens.map(token => token.toLowerCase().trim().replace(/"/g, '')).filter(token => token !== '');
        if (tokens.length > 0) {
            let textToOption = new Map();
            for (let i = 0; i < props.options.length; i++) {
                const option = props.options[i];
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
        props.onChange(event, results);
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
    };

    function showDragIndicator(show) {
        ref.current.style.background = show ? '#1976d2' : '';
    }

    function onDragOver(event) {
        event.preventDefault();
        event.stopPropagation();
        showDragIndicator(true);

    };

    function onDragEnd(event) {
        event.preventDefault();
        event.stopPropagation();
        showDragIndicator(false);
    };


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

    const renderOption = (option, {inputValue}) => {
        inputValue = inputValue.toLowerCase();
        const text = getOptionLabel(option);
        if (inputValue !== '') {
            const index = text.toLowerCase().indexOf(inputValue);
            if (index !== -1) {
                const inputLength = inputValue.length;
                const before = text.substring(0, index);
                const match = text.substring(index, index + inputLength);
                const after = text.substring(index + inputLength);
                return <Typography title={text} noWrap>
                    {before}<b>{match}</b>{after}
                </Typography>;
            }
        }
        return <Typography title={text} noWrap>{text}</Typography>;
    };
    // const renderOption = props.groupBy ? (option) => {
    //  const text = getOptionLabel(option);
    //    return <Typography title={text}  noWrap>{text}</Typography>;
    // };

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
            }
        }
        return exactMatches.concat(startsWithMatches).concat(containsMatches);
    };
    return (
        <Autocomplete
            multiple
            ref={ref}
            filterOptions={filterOptions}
            disableListWrap
            classes={classes}
            getOptionSelected={getOptionSelected}
            value={props.value}
            openOnFocus={true}
            filterSelectedOptions={true}
            getOptionLabel={getOptionLabel}
            groupBy={props.groupBy ? (option) => option.group : null}
            blurOnSelect={true}
            ChipProps={{size: 'small'}}
            ListboxComponent={ListboxComponent}
            renderGroup={renderGroup}
            options={props.options}
            autoHighlight={true}
            onChange={props.onChange}
            renderTags={(value, getTagProps) =>
                value.map((option, index) => {
                    return (
                        <Chip
                            variant="default"
                            onClick={onChipClick ? event => onChipClick(event, option) : null}
                            label={getChipText(option)}
                            title={getChipTitle(option)}
                            size="small"
                            icon={getChipIcon(option)}
                            {...getTagProps({index})}
                        />);
                })
            }
            renderInput={(params) => <TextField margin="dense" {...params} label={props.label}
                                                helperText={props.helperText}/>}
            renderOption={renderOption}
            onPaste={onPaste}
            onDrop={onDrop}
            onDragOver={onDragOver}
            onDragEnd={onDragEnd}
            onDragLeave={onDragEnd}
        />
    );
}

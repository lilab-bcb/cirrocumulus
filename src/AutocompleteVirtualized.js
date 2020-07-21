import {Typography} from '@material-ui/core';
import ListSubheader from '@material-ui/core/ListSubheader';
import {makeStyles, useTheme} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import useMediaQuery from '@material-ui/core/useMediaQuery';
import Autocomplete, {createFilterOptions} from '@material-ui/lab/Autocomplete';
import PropTypes from 'prop-types';
import React from 'react';
import {VariableSizeList} from 'react-window';

const LISTBOX_PADDING = 8; // px

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
    const itemSize = smUp ? 36 : 48;

    const getChildSize = (child) => {
        if (React.isValidElement(child) && child.type === ListSubheader) {
            return 48;
        }

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
                    height={getHeight() + 2 * LISTBOX_PADDING}
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
    <ListSubheader component="div">
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
            let text = event.target.result;
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
        let results = [];
        tokens.forEach(token => {
            token = token.toLowerCase().trim().replace(/"/g, '');
            if (token !== '') {
                for (let i = 0; i < props.options.length; i++) {
                    if (props.options[i].toLowerCase() === token) {
                        results.push(props.options[i]);
                        break;
                    }

                }
            }
        });
        props.onChange(event, results);
    }

    function onPaste(event) {
        let text = event.clipboardData.getData('text/plain');
        if (text != null && text.length > 0) {
            event.preventDefault();
            event.stopPropagation();
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

    const filterOptions = createFilterOptions({matchFrom: 'start'});
    return (
        <Autocomplete
            multiple
            ref={ref}
            filterOptions={filterOptions}
            disableListWrap
            classes={classes}
            getOptionSelected={props.groupBy ? (option, value) => option.text === value.text && option.group === value.group : (option, value) => option === value}
            value={props.value}
            openOnFocus={true}
            filterSelectedOptions={true}
            getOptionLabel={props.groupBy ? (option) => option.text : (option) => option}
            groupBy={props.groupBy ? (option) => option.group : null}
            blurOnSelect={true}
            ChipProps={{size: 'small'}}
            ListboxComponent={ListboxComponent}
            renderGroup={renderGroup}
            options={props.options}
            onChange={props.onChange}
            renderInput={(params) => <TextField {...params} label={props.label}/>}
            renderOption={props.groupBy ? (option) => <Typography title={option.text}
                                                                  noWrap>{option.text}</Typography> : (option) =>
                <Typography
                    noWrap>{option}</Typography>}
            onPaste={onPaste}
            onDrop={onDrop}
            onDragOver={onDragOver}
            onDragEnd={onDragEnd}
            onDragLeave={onDragEnd}
        />
    );
}

import Accordion from '@mui/material/Accordion';
import AccordionDetails from '@mui/material/AccordionDetails';
import AccordionSummary from '@mui/material/AccordionSummary';
import withStyles from '@mui/styles/withStyles';

export const AccordionStyled = withStyles({
    root: {
        border: '1px solid rgba(0, 0, 0, .125)',
        boxShadow: 'none',
        '&:not(:last-child)': {
            borderBottom: 0,
        },
        '&:before': {
            display: 'none',
        },
        '&$expanded': {
            margin: 0,
        },
    },
    expanded: {},
})(Accordion);

export const AccordionSummaryStyled = withStyles({
    root: {
        backgroundColor: 'rgba(0, 0, 0, .03)',
        borderBottom: '1px solid rgba(0, 0, 0, .125)',
        marginBottom: -1,
        minHeight: 43,
        '&$expanded': {
            minHeight: 43,
        },
    },

    content: {
        '&$expanded': {
            margin: 0,
        },
    },
    expanded: {},
})(AccordionSummary);

export const AccordionDetailsStyled = withStyles(theme => ({
    root: {
        flexDirection: 'column',
        padding: 0,
    },
}))(AccordionDetails);

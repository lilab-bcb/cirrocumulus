import {Divider, Paper} from '@mui/material';
import Box from '@mui/material/Box';
import Link from '@mui/material/Link';
import Typography from '@mui/material/Typography';
import ReactMarkdown from 'markdown-to-jsx';
import React from 'react';
import {connect} from 'react-redux';
import {REACT_MD_OVERRIDES} from './util';
import preval from 'preval.macro';

function LandingPage(props) {
    const buildDate = preval`module.exports = new Date().toDateString()`;
    return <Paper elevation={0}>

        <p>Cirrocumulus: an interactive tool for large-scale single-cell genomics</p>
        <Typography variant="h5">Links</Typography>
        <ul>
            <li><Link target="_blank" rel="noopener noreferrer"
                      href="http://cirrocumulus.readthedocs.io/">Documentation</Link></li>
            <li><Link target="_blank" rel="noopener noreferrer"
                      href="https://github.com/klarman-cell-observatory/cirrocumulus">Source
                Code</Link></li>
        </ul>
        <Typography variant="h5">Primary Embedding Controls</Typography>
        <ul>
            <li>Pan: Mouse left click (2-d), right click (3-d)</li>
            <li>Rotate 3-d: Mouse left click</li>
            <li>Zoom: Mouse wheel</li>
            <li>Select: When using lasso or select tool, hold down the Ctrl or Command key to add to selection</li>
            <li>Resize: Click and drag the divider below the primary embedding</li>
            <li>Tooltip: Mouse move</li>
            <li>Select Category: Mouse click</li>
        </ul>

        <Typography variant="h5">Embedding Gallery</Typography>
        <ul>
            <li>Drag charts to reorder</li>
            <li>Click chart to set primary view</li>
        </ul>
        {props.serverInfo && props.serverInfo.footer &&
        <Box><ReactMarkdown options={{overrides: REACT_MD_OVERRIDES}} children={props.serverInfo.footer}/></Box>}

        <Divider/>
        {process.env.REACT_APP_VERSION != null &&
        <Typography variant="body2">Version: {process.env.REACT_APP_VERSION}</Typography>}
        <Typography variant="caption" display="block">{buildDate}</Typography>
    </Paper>;
}


const mapStateToProps = state => {
    return {
        serverInfo: state.serverInfo
    };
};
const mapDispatchToProps = dispatch => {
    return {};
};

export default (connect(
    mapStateToProps, mapDispatchToProps
)(LandingPage));

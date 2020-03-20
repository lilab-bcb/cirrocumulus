import {Tooltip} from '@material-ui/core';
import IconButton from '@material-ui/core/IconButton';
import withStyles from '@material-ui/core/styles/withStyles';
import PauseIcon from '@material-ui/icons/Pause';
import PlayArrowIcon from '@material-ui/icons/PlayArrow';
import React from 'react';

const styles = theme => ({

    root: {
        '& > *': {
            margin: theme.spacing(.4),
        },
        '& > .MuiIconButton-root': {
            padding: 0,
        },
        '& > .cirro-active': {
            fill: 'black'
        },
        '& > .cirro-inactive': {
            fill: '#bdbdbd'
        }
    }
});

//"zoom" | "pan" | "select" | "lasso" | "orbit" | "turntable"
class ChartToolbar extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {dragmode: 'select'};
    }

    setDragMode = (mode) => {
        this.setState({dragmode: mode});
        this.props.onDragMode(mode);
    };

    render() {
        let active = 'cirro-active';
        let inactive = 'cirro-inactive';
        let dragmode = this.state.dragmode;
        return (<div className={this.props.classes.root}>
            {/*<Tooltip title={"Lasso"}>*/}
            {/*    <IconButton className={dragmode === 'lasso' ? active : inactive}*/}
            {/*                aria-label="Lasso Select" onClick={() => this.setDragMode('lasso')}>*/}
            {/*        <svg width="16" height="16" viewBox="0 0 1031 1000">*/}
            {/*            <path*/}
            {/*                d="m1018 538c-36 207-290 336-568 286-277-48-473-256-436-463 10-57 36-108 76-151-13-66 11-137 68-183 34-28 75-41 114-42l-55-70 0 0c-2-1-3-2-4-3-10-14-8-34 5-45 14-11 34-8 45 4 1 1 2 3 2 5l0 0 113 140c16 11 31 24 45 40 4 3 6 7 8 11 48-3 100 0 151 9 278 48 473 255 436 462z m-624-379c-80 14-149 48-197 96 42 42 109 47 156 9 33-26 47-66 41-105z m-187-74c-19 16-33 37-39 60 50-32 109-55 174-68-42-25-95-24-135 8z m360 75c-34-7-69-9-102-8 8 62-16 128-68 170-73 59-175 54-244-5-9 20-16 40-20 61-28 159 121 317 333 354s407-60 434-217c28-159-121-318-333-355z"*/}
            {/*                transform="matrix(1 0 0 -1 0 850)"></path>*/}
            {/*        </svg>*/}
            {/*    </IconButton>*/}
            {/*</Tooltip>*/}
            <Tooltip title={"Brush"}>
                <IconButton edge={false} size={'small'} className={dragmode === 'select' ? active : inactive}
                            aria-label="Box Select" onClick={() => this.setDragMode('select')}>
                    <svg viewBox="0 0 1000 1000" height="16" width="16">
                        <path
                            d="m0 850l0-143 143 0 0 143-143 0z m286 0l0-143 143 0 0 143-143 0z m285 0l0-143 143 0 0 143-143 0z m286 0l0-143 143 0 0 143-143 0z m-857-286l0-143 143 0 0 143-143 0z m857 0l0-143 143 0 0 143-143 0z m-857-285l0-143 143 0 0 143-143 0z m857 0l0-143 143 0 0 143-143 0z m-857-286l0-143 143 0 0 143-143 0z m286 0l0-143 143 0 0 143-143 0z m285 0l0-143 143 0 0 143-143 0z m286 0l0-143 143 0 0 143-143 0z"
                            transform="matrix(1 0 0 -1 0 850)"></path>
                    </svg>
                </IconButton>
            </Tooltip>
            <Tooltip title={"Pan"}>
                <IconButton edge={false} size={'small'} className={dragmode === 'pan' ? active : inactive}
                            aria-label="Pan" onClick={() => this.setDragMode('pan')}>
                    <svg viewBox="0 0 1000 1000" height="16" width="16">
                        <path
                            d="m1000 350l-187 188 0-125-250 0 0 250 125 0-188 187-187-187 125 0 0-250-250 0 0 125-188-188 186-187 0 125 252 0 0-250-125 0 187-188 188 188-125 0 0 250 250 0 0-126 187 188z"
                            transform="matrix(1 0 0 -1 0 850)"></path>
                    </svg>
                </IconButton>
            </Tooltip>
            {this.props.is3d && <Tooltip title={this.props.animating ? 'Pause' : 'Animate'}>
                <IconButton edge={false} size={'small'} aria-label={this.props.animating ? 'Pause' : 'Animate'}
                            onClick={this.props.toggleAnimation}>
                    {!this.props.animating && <PlayArrowIcon/>}
                    {this.props.animating && <PauseIcon/>}
                </IconButton>
            </Tooltip>}


            {/*<Tooltip title={"Zoom"}>*/}
            {/*    <IconButton className={dragmode === 'zoom' ? active : inactive}*/}
            {/*                aria-label="Zoom" onClick={() => this.setDragMode('zoom')}>*/}
            {/*        <svg viewBox="0 0 1000 1000" height="16" width="16">*/}
            {/*            <path*/}
            {/*                d="m1000-25l-250 251c40 63 63 138 63 218 0 224-182 406-407 406-224 0-406-182-406-406s183-406 407-406c80 0 155 22 218 62l250-250 125 125z m-812 250l0 438 437 0 0-438-437 0z m62 375l313 0 0-312-313 0 0 312z"*/}
            {/*                transform="matrix(1 0 0 -1 0 850)"></path>*/}
            {/*        </svg>*/}
            {/*    </IconButton>*/}
            {/*</Tooltip>*/}
            {/*<IconButton className={dragmode === ChartToolbar.MODE_ZOOM_OUT ? active : inactive}*/}
            {/*            aria-label="Zoom Out" onClick={this.onZoomOut}>*/}
            {/*    <svg viewBox="0 0 875 1000" height="16" width="16">*/}
            {/*        <path d="m0 788l0-876 875 0 0 876-875 0z m688-500l-500 0 0 125 500 0 0-125z"*/}
            {/*              transform="matrix(1 0 0 -1 0 850)"></path>*/}
            {/*    </svg>*/}
            {/*</IconButton>*/}

            {!this.props.is3d && <Tooltip title={"Save Image"}>
                <IconButton aria-label="Save Image" onClick={this.props.onSaveImage}>
                    <svg viewBox="0 0 1000 1000" width="16" height="16">
                        <path
                            d="m500 450c-83 0-150-67-150-150 0-83 67-150 150-150 83 0 150 67 150 150 0 83-67 150-150 150z m400 150h-120c-16 0-34 13-39 29l-31 93c-6 15-23 28-40 28h-340c-16 0-34-13-39-28l-31-94c-6-15-23-28-40-28h-120c-55 0-100-45-100-100v-450c0-55 45-100 100-100h800c55 0 100 45 100 100v450c0 55-45 100-100 100z m-400-550c-138 0-250 112-250 250 0 138 112 250 250 250 138 0 250-112 250-250 0-138-112-250-250-250z m365 380c-19 0-35 16-35 35 0 19 16 35 35 35 19 0 35-16 35-35 0-19-16-35-35-35z"
                            transform="matrix(1 0 0 -1 0 850)"></path>
                    </svg>
                </IconButton>
            </Tooltip>}

            {/*<Tooltip title={"Reset"}>*/}
            {/*    <IconButton aria-label="Reset" onClick={this.props.onHome}>*/}
            {/*        <svg viewBox="0 0 928.6 1000" height="16" width="16">*/}
            {/*            <path*/}
            {/*                d="m786 296v-267q0-15-11-26t-25-10h-214v214h-143v-214h-214q-15 0-25 10t-11 26v267q0 1 0 2t0 2l321 264 321-264q1-1 1-4z m124 39l-34-41q-5-5-12-6h-2q-7 0-12 3l-386 322-386-322q-7-4-13-4-7 2-12 7l-35 41q-4 5-3 13t6 12l401 334q18 15 42 15t43-15l136-114v109q0 8 5 13t13 5h107q8 0 13-5t5-13v-227l122-102q5-5 6-12t-4-13z"*/}
            {/*                transform="matrix(1 0 0 -1 0 850)"></path>*/}

            {/*        </svg>*/}
            {/*    </IconButton>*/}
            {/*</Tooltip>*/}
        </div>);
    }
}


export default withStyles(styles)(ChartToolbar);

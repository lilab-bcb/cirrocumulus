import {SvgIcon, Tooltip} from '@mui/material';
import IconButton from '@mui/material/IconButton';
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';

import BorderInnerIcon from '@mui/icons-material/BorderInner';
import CloudQueueIcon from '@mui/icons-material/CloudQueue';
import ExposureIcon from '@mui/icons-material/Exposure';
import FontDownloadIcon from '@mui/icons-material/FontDownload';
import HomeIcon from '@mui/icons-material/Home';
import PauseIcon from '@mui/icons-material/Pause';
import PhotoCameraIcon from '@mui/icons-material/PhotoCamera';
import PhotoLibraryIcon from '@mui/icons-material/PhotoLibrary';
import PhotoSizeSelectSmallIcon from '@mui/icons-material/PhotoSizeSelectSmall';
import PlayArrowIcon from '@mui/icons-material/PlayArrow';
import SettingsIcon from '@mui/icons-material/Settings';
import ZoomInIcon from '@mui/icons-material/ZoomIn';
import ZoomOutIcon from '@mui/icons-material/ZoomOut';
import React, {useState} from 'react';

const active = 'cirro-active';


//"zoom" | "pan" | "select" | "lasso" | "orbit" | "turntable"
export default function ChartToolbar(props) {

    const [saveImageEl, setSaveImageEl] = useState(null);


    function setDragMode(mode) {
        props.onDragMode(mode);
    }

    function onShowLabels() {
        props.onShowLabels();
    }

    function onShowAxis() {
        props.onShowAxis();
    }


    function onShowFog() {
        props.onShowFog();
    }


    function onGallery() {
        props.onGallery();
    }

    function handleMoreOptions() {
        props.onMoreOptions();
    }

    function setEditSelection() {
        props.onEditSelection();
    }

    function handleSaveImageMenu(event) {
        setSaveImageEl(event.currentTarget);

    }

    function handleSaveImageMenuClose(event) {
        setSaveImageEl(null);
    }

    function handleSaveImage(format) {
        setSaveImageEl(null);
        props.onSaveImage(format);
    }


    const {dragmode, editSelection, showLabels, showAxis} = props;
    return <>

        {props.onZoomIn && <Tooltip title={"Zoom In"}>
            <IconButton edge={false} size={'small'}
                        aria-label="Zoom In" onClick={props.onZoomIn}>
                <ZoomInIcon/>
            </IconButton>
        </Tooltip>}

        {props.onZoomOut && <Tooltip title={"Zoom Out"}>
            <IconButton edge={false} size={'small'}
                        aria-label="Zoom In" onClick={props.onZoomOut}>
                <ZoomOutIcon/>
            </IconButton>
        </Tooltip>}

        {props.onHome && <Tooltip title={"Home"}>
            <IconButton edge={false} size={'small'}
                        aria-label="Home" onClick={props.onHome}>
                <HomeIcon/>
            </IconButton>
        </Tooltip>}

        {props.onDragMode && <Tooltip title={"Lasso Select"}>
            <IconButton edge={false} size={'small'} className={dragmode === 'lasso' ? active : ''}
                        aria-label="Lasso Select" onClick={() => setDragMode('lasso')}>
                <SvgIcon width="24" height="21" viewBox="0 0 1031 1000">
                    <path
                        d="m1018 538c-36 207-290 336-568 286-277-48-473-256-436-463 10-57 36-108 76-151-13-66 11-137 68-183 34-28 75-41 114-42l-55-70 0 0c-2-1-3-2-4-3-10-14-8-34 5-45 14-11 34-8 45 4 1 1 2 3 2 5l0 0 113 140c16 11 31 24 45 40 4 3 6 7 8 11 48-3 100 0 151 9 278 48 473 255 436 462z m-624-379c-80 14-149 48-197 96 42 42 109 47 156 9 33-26 47-66 41-105z m-187-74c-19 16-33 37-39 60 50-32 109-55 174-68-42-25-95-24-135 8z m360 75c-34-7-69-9-102-8 8 62-16 128-68 170-73 59-175 54-244-5-9 20-16 40-20 61-28 159 121 317 333 354s407-60 434-217c28-159-121-318-333-355z"
                        transform="matrix(1 0 0 -1 0 850)"></path>

                </SvgIcon>
            </IconButton>
        </Tooltip>}

        {props.onDragMode && <Tooltip title={"Box Select"}>
            <IconButton edge={false} size={'small'}
                        aria-label="Box Select" onClick={() => setDragMode('select')}>
                <PhotoSizeSelectSmallIcon className={dragmode === 'select' ? active : ''}/>
            </IconButton>
        </Tooltip>}

        {props.onEditSelection && <Tooltip title={"Append to selection"}>
            <IconButton edge={false} size={'small'} aria-label="Append to selection"
                        onClick={setEditSelection}>
                <ExposureIcon className={editSelection ? active : ''}/>
            </IconButton>
        </Tooltip>}

        {props.is3d && <Tooltip title={props.animating ? 'Pause' : 'Animate'}>
            <IconButton edge={false} size={'small'}
                        aria-label={props.animating ? 'Pause' : 'Animate'}
                        onClick={props.toggleAnimation}>
                {!props.animating && <PlayArrowIcon/>}
                {props.animating && <PauseIcon/>}
            </IconButton>
        </Tooltip>}

        {props.onDragMode && <Tooltip title={"Pan"}>
            <IconButton edge={false} size={'small'} className={dragmode === 'pan' ? active : ''}
                        aria-label="Pan" onClick={() => setDragMode('pan')}>
                <SvgIcon viewBox="0 0 1000 1000" height="16" width="16">
                    <path
                        d="m1000 350l-187 188 0-125-250 0 0 250 125 0-188 187-187-187 125 0 0-250-250 0 0 125-188-188 186-187 0 125 252 0 0-250-125 0 187-188 188 188-125 0 0 250 250 0 0-126 187 188z"
                        transform="matrix(1 0 0 -1 0 850)"></path>
                </SvgIcon>
            </IconButton>
        </Tooltip>}

        {props.onShowLabels && <Tooltip title={"Show Categorical Labels"}>
            <IconButton edge={false} size={'small'} className={showLabels ? active : ''}
                        aria-label="Show Labels" onClick={() => onShowLabels()}>
                <FontDownloadIcon/>
            </IconButton>
        </Tooltip>}


        {props.is3d && props.onShowAxis && <Tooltip title={"Show Axis"}>
            <IconButton edge={false} size={'small'} className={showAxis ? active : ''}
                        aria-label="Show Axis" onClick={() => onShowAxis()}>
                <BorderInnerIcon/>
            </IconButton>
        </Tooltip>}


        {props.is3d && props.onShowFog && <Tooltip title={"Show Fog"}>
            <IconButton edge={false} size={'small'} className={props.showFog ? active : ''}
                        aria-label="Show Fog" onClick={() => onShowFog()}>
                <CloudQueueIcon/>
            </IconButton>
        </Tooltip>}


        {/*<Tooltip title={"Zoom"}>*/}
        {/*    <IconButton className={dragmode === 'zoom' ? active : ''}*/}
        {/*                aria-label="Zoom" onClick={() => setDragMode('zoom')}>*/}
        {/*        <svg viewBox="0 0 1000 1000" height="16" width="16">*/}
        {/*            <path*/}
        {/*                d="m1000-25l-250 251c40 63 63 138 63 218 0 224-182 406-407 406-224 0-406-182-406-406s183-406 407-406c80 0 155 22 218 62l250-250 125 125z m-812 250l0 438 437 0 0-438-437 0z m62 375l313 0 0-312-313 0 0 312z"*/}
        {/*                transform="matrix(1 0 0 -1 0 850)"></path>*/}
        {/*        </svg>*/}
        {/*    </IconButton>*/}
        {/*</Tooltip>*/}
        {/*<IconButton className={dragmode === ChartToolbar.MODE_ZOOM_OUT ? active : ''}*/}
        {/*            aria-label="Zoom Out" onClick={onZoomOut}>*/}
        {/*    <svg viewBox="0 0 875 1000" height="16" width="16">*/}
        {/*        <path d="m0 788l0-876 875 0 0 876-875 0z m688-500l-500 0 0 125 500 0 0-125z"*/}
        {/*              transform="matrix(1 0 0 -1 0 850)"></path>*/}
        {/*    </svg>*/}
        {/*</IconButton>*/}

        <Tooltip title={"Save Image"}>
            <IconButton aria-controls="save-image-menu" aria-haspopup="true" edge={false} size={'small'}
                        aria-label="Save Image" onClick={handleSaveImageMenu}>
                <PhotoCameraIcon/>
            </IconButton>
        </Tooltip>

        <Menu
            id="save-image-menu"
            anchorEl={saveImageEl}
            keepMounted
            open={Boolean(saveImageEl)}
            onClose={handleSaveImageMenuClose}
        >
            <MenuItem onClick={e => handleSaveImage('png')}>PNG</MenuItem>
            <MenuItem onClick={e => handleSaveImage('svg')}>SVG</MenuItem>

        </Menu>


        {props.onMoreOptions && <Tooltip title={"More Options"}>
            <IconButton edge={false} size={'small'}
                        aria-label="More Options" onClick={handleMoreOptions}>
                <SettingsIcon/>
            </IconButton>
        </Tooltip>}


        <Tooltip title={"Scroll To Gallery"}>
            <IconButton edge={false} size={'small'}
                        aria-label="Scroll To Gallery" onClick={onGallery}>
                <PhotoLibraryIcon/>
            </IconButton>
        </Tooltip>
        {/*<Tooltip title={"Copy Image"}>*/}
        {/*    <IconButton edge={false} size={'small'} aria-label="Copy Image" onClick={props.onCopyImage}>*/}
        {/*        <FileCopyIcon/>*/}
        {/*    </IconButton>*/}
        {/*</Tooltip>*/}

        {/*<Tooltip title={*/}
        {/*    <>*/}
        {/*        <h6>3D controls</h6>*/}
        {/*        <b>Rotate</b> Mouse left click.*/}
        {/*        <b>Pan</b> Mouse right click.*/}
        {/*        <b>Zoom</b> Mouse wheel.*/}
        {/*        Holding <b>ctrl</b> reverses the mouse clicks.*/}
        {/*        <h6>2D controls</h6>*/}
        {/*        <b>Pan</b> Mouse left click.*/}
        {/*        <b>Zoom</b> Mouse wheel.*/}
        {/*    </>*/}
        {/*}>*/}
        {/*    <IconButton edge={false} size={'small'} aria-label="Help">*/}
        {/*        <HelpIcon/>*/}
        {/*    </IconButton>*/}
        {/*</Tooltip>*/}
        {/*<Tooltip title={"Reset"}>*/}
        {/*    <IconButton aria-label="Reset" onClick={props.onHome}>*/}
        {/*        <svg viewBox="0 0 928.6 1000" height="16" width="16">*/}
        {/*            <path*/}
        {/*                d="m786 296v-267q0-15-11-26t-25-10h-214v214h-143v-214h-214q-15 0-25 10t-11 26v267q0 1 0 2t0 2l321 264 321-264q1-1 1-4z m124 39l-34-41q-5-5-12-6h-2q-7 0-12 3l-386 322-386-322q-7-4-13-4-7 2-12 7l-35 41q-4 5-3 13t6 12l401 334q18 15 42 15t43-15l136-114v109q0 8 5 13t13 5h107q8 0 13-5t5-13v-227l122-102q5-5 6-12t-4-13z"*/}
        {/*                transform="matrix(1 0 0 -1 0 850)"></path>*/}

        {/*        </svg>*/}
        {/*    </IconButton>*/}
        {/*</Tooltip>*/}
    </>;

}



/*
@license
Copyright 2019 Google LLC. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/


import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls';

import {CameraType, LabelRenderParams, RenderContext} from './render';
import {Styles} from './styles';
import {InteractionMode, Point2D, Point3D} from './types';
import * as util from './util';

import {ScatterPlotVisualizer} from './scatter_plot_visualizer';
import {Point, ScatterBoundingBox, ScatterPlotRectangleSelector} from './scatter_plot_rectangle_selector';
import {
    AxesHelper,
    Camera,
    MOUSE,
    Object3D,
    OrthographicCamera,
    PerspectiveCamera,
    PointLight,
    Scene,
    Vector3,
    WebGLRenderer,
} from 'three';

/**
 * The length of the cube (diameter of the circumscribing sphere) where all the
 * points live.
 */
const CUBE_LENGTH = 2;
const MAX_ZOOM = 5 * CUBE_LENGTH;
const MIN_ZOOM = 0.025 * CUBE_LENGTH;

// Constants relating to the camera parameters.
const PERSP_CAMERA_FOV_VERTICAL = 70;
const PERSP_CAMERA_NEAR_CLIP_PLANE = 0.01;
const PERSP_CAMERA_FAR_CLIP_PLANE = 100;
const ORTHO_CAMERA_FRUSTUM_HALF_EXTENT = 1.2;

const START_CAMERA_POS_3D = new Vector3(0.45, 0.9, 1.6);
const START_CAMERA_TARGET_3D = new Vector3(0, 0, 0);
const START_CAMERA_POS_2D = new Vector3(0, 0, 4);
const START_CAMERA_TARGET_2D = new Vector3(0, 0, 0);

const ORBIT_MOUSE_ROTATION_SPEED = 1;
const ORBIT_ANIMATION_ROTATION_CYCLE_IN_SECONDS = 7;
const ORBIT_ZOOM_SPEED = 1.125;

export type OnCameraMoveListener = (
    cameraPosition: Vector3,
    cameraTarget: Vector3
) => void;

/** Defines a camera, suitable for serialization. */
export interface CameraDef {
    orthographic: boolean;
    position: Point3D;
    target: Point3D;
    zoom: number;
}

export type CameraParams = Partial<Pick<CameraDef, 'position' | 'target' | 'zoom'>>;

export interface ScatterPlotParams {
    camera?: CameraParams;
    onClick?: (point: number | null) => void;
    onHover?: (point: number | null) => void;
    onSelect?: (points: number[], boundingBox?: ScatterBoundingBox) => void;
    selectEnabled?: boolean;
    styles: Styles;
    interactive?: boolean;
}

/**
 * Maintains a js instantiation and context,
 * animation state, and all other logic that's
 * independent of how a 3D scatter plot is actually rendered. Also holds an
 * array of visualizers and dispatches application events to them.
 */
export class ScatterPlot {
    private container: HTMLElement;
    private styles: Styles;
    private interactive: boolean = true;
    private hoverPoint: Point = {x: 0, y: 0};
    private hoverCallback: (point: Point | null, e: PointerEvent | null) => void = () => {
    };
    private clickCallback: (point: Point | null, appendToSelection: boolean) => void = () => {
    };
    private boxCallback: (boundingBox: ScatterBoundingBox, appendToSelection: boolean) => void = () => {
    };
    private lassoCallback: (points: Point[], appendToSelection: boolean) => void = () => {
    };

    private cameraCallback: (eventName: string, cameraPosition: Vector3, cameraTarget: Vector3) => void = () => {
    };
    // Array of visualizers
    private visualizers: ScatterPlotVisualizer[] = [];

    private height = 0;
    private width = 0;
    private dimensions = 3;

    private interactionMode = InteractionMode.PAN;

    private renderer: WebGLRenderer;

    private scene: Scene;

    private light: PointLight;

    private camera!: Camera;
    private orbitAnimationOnNextCameraCreation: boolean = false;
    private orbitCameraControls: any;
    private orbitAnimationId: number | null = null;

    private worldSpacePointPositions = new Float32Array(0);
    private pointColors = new Float32Array(0);
    private pointScaleFactors = new Float32Array(0);
    private labels?: LabelRenderParams;
    private polylineColors: { [polylineIndex: number]: Float32Array } = {};
    private polylineOpacities = new Float32Array(0);
    private polylineWidths = new Float32Array(0);

    private selecting = false;
    private mouseIsDown = false;
    private isDragSequence = false;
    private rectangleSelector!: ScatterPlotRectangleSelector;

    constructor(containerElement: HTMLElement, params: ScatterPlotParams, premultipliedAlpha: boolean = false) {
        this.container = containerElement;
        this.styles = params.styles;
        this.interactive = params.interactive || false;

        this.computeLayoutValues();

        this.scene = new Scene();
        this.renderer = new WebGLRenderer({
            alpha: true,
            premultipliedAlpha: premultipliedAlpha,
            antialias: false,
        });


        this.renderer.setClearColor(this.styles.backgroundColor, 1);
        this.container.appendChild(this.renderer.domElement);
        this.light = new PointLight(0xffecbf, 1, 0);
        this.scene.add(this.light);

        if (params.interactive) {
            this.rectangleSelector = new ScatterPlotRectangleSelector(
                this.container,
                this.selectBoundingBox.bind(this),
                this.selectLasso.bind(this),
                this.styles
            );
            this.addInteractionListeners();
        }
        this.setDimensions(3);
        this.makeCamera(params.camera);
        this.resize();
    }

    private addInteractionListeners() {
        this.container.addEventListener('pointerout', this.onPointerOut.bind(this));
        this.container.addEventListener('pointerenter', this.onPointerOut.bind(this));
        this.container.addEventListener('pointermove', this.onPointerMove.bind(this));
        this.container.addEventListener('pointerdown', this.onPointerDown.bind(this));
        this.container.addEventListener('pointerup', this.onPointerUp.bind(this));
        // this.container.addEventListener('click', this.onClick.bind(this));
        window.addEventListener('keydown', this.onKeyDown.bind(this), false);
        window.addEventListener('keyup', this.onKeyUp.bind(this), false);
    }

    private addCameraControlsEventListeners(cameraControls: any) {
        if (!this.interactive) {
            return;
        }
        // Start is called when the user stars interacting with
        // controls.
        cameraControls.addEventListener('start', () => {
            this.stopOrbitAnimation();
            this.cameraCallback('start', this.camera.position, cameraControls.target);
        });

        // Change is called everytime the user interacts with the controls.
        cameraControls.addEventListener('change', () => {
            this.render();
            this.cameraCallback('change', this.camera.position, cameraControls.target);
        });

        // End is called when the user stops interacting with the
        // controls (e.g. on mouse up, after dragging).
        cameraControls.addEventListener('end', () => {
            this.cameraCallback('end', this.camera.position, cameraControls.target);
        });
    }

    private makeOrbitControls(camera: Camera, cameraIs3D: boolean) {
        if (this.orbitCameraControls != null) {
            this.orbitCameraControls.dispose();
        }

        const occ = new OrbitControls(camera, this.renderer.domElement);

        occ.zoomSpeed = ORBIT_ZOOM_SPEED;
        occ.enableRotate = cameraIs3D;
        occ.autoRotate = false;
        occ.enableKeys = false;
        occ.rotateSpeed = ORBIT_MOUSE_ROTATION_SPEED;
        if (cameraIs3D) {
            occ.mouseButtons.LEFT = MOUSE.LEFT; // Orbit
            occ.mouseButtons.RIGHT = MOUSE.RIGHT; // Pan
        } else {
            occ.mouseButtons.LEFT = MOUSE.RIGHT; // Orbit
            occ.mouseButtons.RIGHT = MOUSE.LEFT; //Pan
        }
        occ.reset();

        this.camera = camera;
        this.orbitCameraControls = occ;

        this.addCameraControlsEventListeners(this.orbitCameraControls);
    }

    private makeCamera(cameraParams: CameraParams = {}) {
        const def = this.makeDefaultCameraDef(this.dimensions, cameraParams);
        this.recreateCamera(def);

        if (this.dimensions === 3 && this.styles.axesVisible) {
            this.add3dAxes();
        } else {
            this.remove3dAxesFromScene();
        }
    }

    private makeCamera3D(cameraDef: CameraDef, w: number, h: number) {
        let camera: PerspectiveCamera;
        {
            const aspectRatio = w / h;
            camera = new PerspectiveCamera(
                PERSP_CAMERA_FOV_VERTICAL,
                aspectRatio,
                PERSP_CAMERA_NEAR_CLIP_PLANE,
                PERSP_CAMERA_FAR_CLIP_PLANE
            );
            camera.position.set(
                cameraDef.position[0],
                cameraDef.position[1],
                cameraDef.position[2]
            );
            const at = new Vector3(
                cameraDef.target[0],
                cameraDef.target[1],
                cameraDef.target[2]
            );
            camera.lookAt(at);
            camera.zoom = cameraDef.zoom;
            camera.updateProjectionMatrix();
        }
        this.camera = camera;

        this.makeOrbitControls(camera, true);

    }

    private makeCamera2D(cameraDef: CameraDef, w: number, h: number) {
        let camera: OrthographicCamera;
        const target = new Vector3(
            cameraDef.target[0],
            cameraDef.target[1],
            cameraDef.target[2]
        );
        {
            const aspectRatio = w / h;
            let left = -ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
            let right = ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
            let bottom = -ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
            let top = ORTHO_CAMERA_FRUSTUM_HALF_EXTENT;
            // Scale up the larger of (w, h) to match the aspect ratio.
            if (aspectRatio > 1) {
                left *= aspectRatio;
                right *= aspectRatio;
            } else {
                top /= aspectRatio;
                bottom /= aspectRatio;
            }
            camera = new OrthographicCamera(
                left,
                right,
                top,
                bottom,
                -1000,
                1000
            );
            camera.position.set(
                cameraDef.position[0],
                cameraDef.position[1],
                cameraDef.position[2]
            );
            // The orbit controls pan up operation is tied to the z dimension
            camera.up = new Vector3(0, 0, 1);
            camera.lookAt(target);
            camera.zoom = cameraDef.zoom;
            camera.updateProjectionMatrix();
        }
        this.camera = camera;
        this.makeOrbitControls(camera, false);
    }

    private makeDefaultCameraDef(
        dimensions: number,
        cameraParams: CameraParams = {}
    ): CameraDef {
        const orthographic = dimensions === 2;
        const position = orthographic ? START_CAMERA_POS_2D : START_CAMERA_POS_3D;
        const target = orthographic
            ? START_CAMERA_TARGET_2D
            : START_CAMERA_TARGET_3D;
        const def: CameraDef = {
            orthographic,
            zoom: 1.0,
            position: [position.x, position.y, position.z],
            target: [target.x, target.y, target.z],
        };

        if (cameraParams.zoom) def.zoom = cameraParams.zoom;
        if (cameraParams.position) def.position = cameraParams.position;
        if (cameraParams.target) def.target = cameraParams.target;

        return def;
    }

    /** Recreate the scatter plot camera from a definition structure. */
    recreateCamera(cameraDef: CameraDef) {
        if (cameraDef.orthographic) {
            this.makeCamera2D(cameraDef, this.width, this.height);
        } else {
            this.makeCamera3D(cameraDef, this.width, this.height);
        }
        if (this.interactive) {
            this.orbitCameraControls.minDistance = MIN_ZOOM;
            this.orbitCameraControls.maxDistance = MAX_ZOOM;
            this.orbitCameraControls.update();
            if (this.orbitAnimationOnNextCameraCreation) {
                this.startOrbitAnimation();
            }
        }
    }

    setInteractionMode(interactionMode: InteractionMode) {
        this.interactionMode = interactionMode;
        if (interactionMode === InteractionMode.SELECT) {
            this.selecting = true;
            this.container.style.cursor = 'crosshair';
        } else {
            this.selecting = false;
            this.container.style.cursor = 'default';
        }
        this.orbitCameraControls.enabled = !this.selecting;
    }


    private onPointerDown(e: PointerEvent) {
        this.hoverCallback(null, null);
        this.isDragSequence = false;
        this.mouseIsDown = true;
        if (this.selecting) {
            this.orbitCameraControls.enabled = false;
            this.rectangleSelector.onMouseDown(e.offsetX, e.offsetY, e.ctrlKey || e.metaKey);
            // this.setNearestPointToMouse(e);
        } else if (
            !e.ctrlKey &&
            this.sceneIs3D() &&
            this.orbitCameraControls.mouseButtons.ORBIT === MOUSE.RIGHT
        ) {
            // The user happened to press the ctrl key when the tab was active,
            // unpressed the ctrl when the tab was inactive, and now he/she
            // is back to the projector tab.
            this.orbitCameraControls.mouseButtons.ORBIT = MOUSE.LEFT;
            this.orbitCameraControls.mouseButtons.PAN = MOUSE.RIGHT;
        } else if (
            e.ctrlKey &&
            this.sceneIs3D() &&
            this.orbitCameraControls.mouseButtons.ORBIT === MOUSE.LEFT
        ) {
            // Similarly to the situation above.
            this.orbitCameraControls.mouseButtons.ORBIT = MOUSE.RIGHT;
            this.orbitCameraControls.mouseButtons.PAN = MOUSE.LEFT;
        }
    }

    /** When we stop dragging/zooming, return to normal behavior. */
    private onPointerUp(e: PointerEvent) {
        if (this.selecting) {
            this.rectangleSelector.onMouseUp();
        } else {
            this.hoverPoint.x = e.offsetX;
            this.hoverPoint.y = e.offsetY;
            this.clickCallback(this.hoverPoint, e.metaKey || e.ctrlKey);
        }
        this.mouseIsDown = false;
    }


    /**
     * When the mouse moves, find the nearest point (if any) and send it to the
     * hoverlisteners (usually called from embedding.ts)
     */
    private onPointerMove(e: PointerEvent) {
        this.isDragSequence = this.mouseIsDown;
        // Depending if we're selecting or just navigating, handle accordingly.
        if (this.selecting && this.mouseIsDown) {
            this.rectangleSelector.onMouseMove(e.offsetX, e.offsetY);
            // this.render();
        } else if (!this.mouseIsDown) {
            // this.setNearestPointToMouse(e);
            this.hoverPoint.x = e.offsetX;
            this.hoverPoint.y = e.offsetY;
            this.hoverCallback(this.hoverPoint, e);
        }
    }

    private onPointerOut(e: PointerEvent) {
        if (!this.selecting) {
            this.hoverCallback(null, null);
        }
    }


    /** For using ctrl + left click as right click, and for circle select */
    private onKeyDown(e: KeyboardEvent) {
        // If ctrl is pressed, use left click to orbit
        if (this.sceneIs3D() && e.shiftKey) {
            this.orbitCameraControls.mouseButtons.ORBIT = MOUSE.RIGHT;
            this.orbitCameraControls.mouseButtons.PAN = MOUSE.LEFT;
        }

        // If shift is pressed, start selecting
        // if (e.keyCode === SHIFT_KEY && this.selectEnabled) {
        //     this.selecting = true;
        //     this.container.style.cursor = 'crosshair';
        // }
    }

    /** For using ctrl + left click as right click, and for circle select */
    private onKeyUp(e: KeyboardEvent) {
        if (this.sceneIs3D() && e.shiftKey) {
            this.orbitCameraControls.mouseButtons.ORBIT = MOUSE.LEFT;
            this.orbitCameraControls.mouseButtons.PAN = MOUSE.RIGHT;
        }

        // If shift is released, stop selecting
        // if (e.keyCode === SHIFT_KEY && this.selectEnabled) {
        //     this.selecting = false;
        //     this.container.style.cursor = 'default';
        //     this.render();
        // }
    }


    private selectBoundingBox(boundingBox: ScatterBoundingBox, appendToSelection: boolean) {
        this.boxCallback(boundingBox, appendToSelection);
    }

    private selectLasso(points: Point[], appendToSelection: boolean) {
        this.lassoCallback(points, appendToSelection);
    }


    private computeLayoutValues(): Point2D {
        this.width = this.container.offsetWidth;
        this.height = Math.max(1, this.container.offsetHeight);
        return [this.width, this.height];
    }

    private sceneIs3D(): boolean {
        return this.dimensions === 3;
    }

    private remove3dAxesFromScene(): Object3D | undefined {
        const axes = this.scene.getObjectByName('axes');
        if (axes != null) {
            this.scene.remove(axes);
        }
        return axes;
    }

    private add3dAxes() {
        const axes = new AxesHelper();
        axes.name = 'axes';
        this.scene.add(axes);
    }

    /** Set 2d vs 3d mode. */
    setDimensions(dimensions: number) {
        if (dimensions !== 2 && dimensions !== 3) {
            throw new RangeError('dimensions must be 2 or 3');
        }

        if (this.dimensions !== dimensions) {
            this.dimensions = dimensions;
            this.makeCamera();
        }
    }

    /** Gets the current camera position. */
    getCameraPosition(): Point3D {
        const currPos = this.camera.position;
        return [currPos.x, currPos.y, currPos.z];
    }

    /** Gets the current camera target. */
    getCameraTarget(): Point3D {
        let currTarget = this.orbitCameraControls.target;
        return [currTarget.x, currTarget.y, currTarget.z];
    }

    /** Sets up the camera from given position and target coordinates. */
    setCameraPositionAndTarget(position: Point3D, target: Point3D) {
        this.stopOrbitAnimation();
        this.camera.position.set(position[0], position[1], position[2]);
        this.orbitCameraControls.target.set(target[0], target[1], target[2]);
        this.orbitCameraControls.update();
        this.render();
    }

    /** Starts orbiting the camera around its current lookat target. */
    startOrbitAnimation() {
        if (!this.sceneIs3D()) {
            return;
        }
        if (this.orbitAnimationId != null) {
            this.stopOrbitAnimation();
        }
        this.orbitCameraControls.autoRotate = true;
        this.orbitCameraControls.rotateSpeed = ORBIT_ANIMATION_ROTATION_CYCLE_IN_SECONDS;
        this.updateOrbitAnimation();
    }

    orbitIsAnimating() {
        return this.orbitAnimationId != null;
    }

    private updateOrbitAnimation() {
        this.orbitCameraControls.update();
        this.orbitAnimationId = requestAnimationFrame(() =>
            this.updateOrbitAnimation()
        );
    }

    /** Stops the orbiting animation on the camera. */
    stopOrbitAnimation() {
        this.orbitCameraControls.autoRotate = false;
        this.orbitCameraControls.rotateSpeed = ORBIT_MOUSE_ROTATION_SPEED;
        if (this.orbitAnimationId != null) {
            cancelAnimationFrame(this.orbitAnimationId);
            this.orbitAnimationId = null;
        }
    }

    getActiveVisualizers(): ScatterPlotVisualizer[] {
        return this.visualizers;
    }

    setActiveVisualizers(visualizers: ScatterPlotVisualizer[]) {
        this.visualizers.forEach(visualizer => {
            if (visualizers.indexOf(visualizer) === -1) {
                visualizer.dispose();
            }
        });
        this.visualizers = visualizers;
        this.visualizers.forEach(visualizer => {
            visualizer.setScene(this.scene);
            visualizer.onResize(this.width, this.height);
            if (this.worldSpacePointPositions) {
                visualizer.onPointPositionsChanged(this.worldSpacePointPositions);
            }
        });
    }

    /** Disposes all visualizers attached to this scatter plot. */
    dispose() {
        this.visualizers.forEach(v => v.dispose());
        this.visualizers = [];
    }

    /** Update scatter plot with a new array of packed xyz point positions. */
    setPointPositions(worldSpacePointPositions: Float32Array) {
        this.worldSpacePointPositions = worldSpacePointPositions;
        this.visualizers.forEach(v =>
            v.onPointPositionsChanged(worldSpacePointPositions)
        );
    }

    render() {
        {
            const lightPos = this.camera.position.clone();
            lightPos.x += 1;
            lightPos.y += 1;
            this.light.position.set(lightPos.x, lightPos.y, lightPos.z);
        }

        const cameraType =
            this.camera instanceof PerspectiveCamera
                ? CameraType.Perspective
                : CameraType.Orthographic;

        let cameraSpacePointExtents: [number, number] = [0, 0];
        if (this.worldSpacePointPositions != null) {
            cameraSpacePointExtents = util.getNearFarPoints(
                this.worldSpacePointPositions,
                this.camera.position,
                this.orbitCameraControls.target
            );
        }

        const rc = new RenderContext(
            this.camera,
            cameraType,
            this.orbitCameraControls.target,
            this.width,
            this.height,
            cameraSpacePointExtents[0],
            cameraSpacePointExtents[1],
            this.styles.backgroundColor,
            this.pointColors,
            this.pointScaleFactors,
            this.labels,
            this.polylineColors,
            this.polylineOpacities,
            this.polylineWidths
        );


        // Render second pass to color buffer, to be displayed on the canvas.
        this.visualizers.forEach(v => v.onRender(rc));

        this.renderer.setRenderTarget(null);
        this.renderer.render(this.scene, this.camera);
    }

    /** Set the colors for every data point. (RGB triplets) */
    setPointColors(colors: Float32Array) {
        this.pointColors = colors;
    }

    /** Set the scale factors for every data point. (scalars) */
    setPointScaleFactors(scaleFactors: Float32Array) {
        this.pointScaleFactors = scaleFactors;
    }

    /** Set the labels to rendered */
    setLabels(labels: LabelRenderParams) {
        this.labels = labels;
    }

    /** Set the colors for every data polyline. (RGB triplets) */
    setPolylineColors(colors: { [polylineIndex: number]: Float32Array }) {
        this.polylineColors = colors;
    }

    setPolylineOpacities(opacities: Float32Array) {
        this.polylineOpacities = opacities;
    }

    setPolylineWidths(widths: Float32Array) {
        this.polylineWidths = widths;
    }

    resetZoom() {
        this.recreateCamera(this.makeDefaultCameraDef(this.dimensions));
        this.render();
    }

    setDayNightMode(isNight: boolean) {
        const canvases = this.container.querySelectorAll('canvas');
        const filterValue = isNight ? 'invert(100%)' : '';
        for (let i = 0; i < canvases.length; i++) {
            canvases[i].style.filter = filterValue;
        }
    }

    resize(render = true) {
        const [oldW, oldH] = [this.width, this.height];
        const [newW, newH] = this.computeLayoutValues();

        if (this.dimensions === 3) {
            const camera = this.camera as PerspectiveCamera;
            camera.aspect = newW / newH;
            camera.updateProjectionMatrix();
        } else {
            const camera = this.camera as OrthographicCamera;
            // Scale the ortho frustum by however much the window changed.
            const scaleW = newW / oldW;
            const scaleH = newH / oldH;
            const newCamHalfWidth = ((camera.right - camera.left) * scaleW) / 2;
            const newCamHalfHeight = ((camera.top - camera.bottom) * scaleH) / 2;
            camera.top = newCamHalfHeight;
            camera.bottom = -newCamHalfHeight;
            camera.left = -newCamHalfWidth;
            camera.right = newCamHalfWidth;
            camera.updateProjectionMatrix();
        }

        // Accouting for retina displays.
        const dpr = window.devicePixelRatio || 1;
        this.renderer.setPixelRatio(dpr);
        this.renderer.setSize(newW, newH);


        this.visualizers.forEach(v => v.onResize(newW, newH));

        if (render) {
            this.render();
        }
    }

    updateFromCameraDef(cameraDef: any) {
        if (cameraDef.position) {
            this.camera.position.set(
                cameraDef.position[0],
                cameraDef.position[1],
                cameraDef.position[2]
            );
        }

        if (cameraDef.target) {
            const at = new Vector3(
                cameraDef.target[0],
                cameraDef.target[1],
                cameraDef.target[2]
            );
            this.camera.lookAt(at);
        }


        if (cameraDef.zoom != null) {
            // @ts-ignore
            this.camera.zoom = cameraDef.zoom;
        }
        // @ts-ignore
        this.camera.updateProjectionMatrix();
    }

    getCameraDef() {
        const def: any = {};
        const pos = this.camera.position;
        const tgt = this.orbitCameraControls.target;
        def.position = [pos.x, pos.y, pos.z];
        def.target = [tgt.x, tgt.y, tgt.z];
        def.zoom = (this.camera as any).zoom;
        return def;
    }

}

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


import {ScatterPlotVisualizer} from './scatter_plot_visualizer';
import {RenderContext} from './render';
import {Styles} from "./styles";
import {Scene, Vector3} from "three";

/**
 * Creates and maintains a 2d canvas on top of the GL canvas. All labels, when
 * active, are rendered to the 2d canvas as part of the visible render pass.
 */
export class ScatterPlotVisualizerCanvasLabels
    implements ScatterPlotVisualizer {
    public id = 'CANVAS_LABELS';

    private worldSpacePointPositions = new Float32Array(0);
    private gc: CanvasRenderingContext2D;
    private canvas: HTMLCanvasElement;
    private labelsActive: boolean = true;
    private labelStrings: string[] = [];
    private fillStyle: string = 'black';
    private shadowColor: string = 'white';
    private font: string = 'bold 14px Roboto Condensed';


    constructor(container: HTMLElement, private styles: Styles) {
        this.canvas = document.createElement('canvas');
        container.appendChild(this.canvas);
        this.gc = this.canvas.getContext('2d')!;
        this.canvas.style.position = 'absolute';
        this.canvas.style.left = '0';
        this.canvas.style.top = '0';
        this.canvas.className = 'label';
        this.canvas.style.pointerEvents = 'none';
    }

    private clear() {
        this.gc.clearRect(0, 0, this.canvas.width * window.devicePixelRatio, this.canvas.height * window.devicePixelRatio);
    }

    setLabels(labelStrings: string[], positions: Float32Array) {
        this.labelStrings = labelStrings;
        this.worldSpacePointPositions = positions;
    }

    private draw(rc: RenderContext) {

        const camera = rc.camera;
        const positions = this.worldSpacePointPositions;
        const pos = new Vector3();
        const labelStrings = this.labelStrings;
        const dpr = window.devicePixelRatio;
        const context = this.gc;
        context.scale(dpr, dpr);
        context.fillStyle = this.fillStyle; // this.styles.label3D.color;
        context.font = this.font;
        context.textAlign = "center";
        context.shadowColor = this.shadowColor;
        context.shadowBlur = 4;
        const widthHalf = (this.canvas.width / dpr) / 2;
        const heightHalf = (this.canvas.height / dpr) / 2;
        for (let i = 0, k = 0; i < labelStrings.length; i++, k += 3) {
            pos.x = positions[k];
            pos.y = positions[k + 1];
            pos.z = positions[k + 2];
            pos.project(camera);
            pos.x = (pos.x * widthHalf) + widthHalf;
            pos.y = -(pos.y * heightHalf) + heightHalf;
            context.fillText(labelStrings[i], pos.x, pos.y);
        }
        context.setTransform(1, 0, 0, 1, 0, 0);
    }

    onResize(newWidth: number, newHeight: number) {
        let dpr = window.devicePixelRatio;
        this.canvas.width = newWidth * dpr;
        this.canvas.height = newHeight * dpr;
        this.canvas.style.width = newWidth + 'px';
        this.canvas.style.height = newHeight + 'px';
        this.gc = this.canvas.getContext('2d')!;
    }

    onRender(rc: RenderContext) {
        this.clear();
        if (this.labelsActive) {
            this.draw(rc);
        }

    }

    dispose() {
        // this.clear();
    }

    onPointPositionsChanged(newPositions: Float32Array) { // ignore

    }

    setScene(scene: Scene) {
    }

    onPickingRender(renderContext: RenderContext) {
    }


}

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


import {
    BufferAttribute,
    BufferGeometry,
    Color,
    Fog,
    NormalBlending,
    PerspectiveCamera,
    Points,
    Scene,
    ShaderChunk,
    ShaderMaterial,
    Texture
} from 'three';
import {ScatterPlotVisualizer} from './scatter_plot_visualizer';
import {RenderContext} from './render';
import {Styles} from './styles';
import * as util from './util';
import {INDEX_NUM_ELEMENTS, RGBA_NUM_ELEMENTS, XYZ_NUM_ELEMENTS,} from './constants';

export interface SpriteSheetParams {
    spritesheetImage: HTMLImageElement | string;
    spriteDimensions: [number, number];
    spriteIndices: Float32Array;
    onImageLoad: () => void;
}

const VERTEX_SHADER = `
    // Index of the specific vertex (passed in as bufferAttribute), and the
    // variable that will be used to pass it to the fragment shader.
    attribute float spriteIndex;
    attribute vec4 color;
    attribute float scaleFactor;

    varying vec4 vColor;

    uniform bool sizeAttenuation;
    uniform float pointSize;
    uniform float zoomFactor;
    uniform bool zoomFactorSpecified;
    
    varying float fogDepth;

    void main() {
      // Pass index and color values to fragment shader.
      vColor = color;

      // Transform current vertex by modelViewMatrix (model world position and
      // camera world position matrix).
      vec4 cameraSpacePos = modelViewMatrix * vec4(position, 1.0);

      // Project vertex in camera-space to screen coordinates using the camera's
      // projection matrix.
      gl_Position = projectionMatrix * cameraSpacePos;

      // Create size attenuation (if we're in 3D mode) by making the size of
      // each point inversely proportional to its distance to the camera.
      float outputPointSize = pointSize;
      if (sizeAttenuation) {
        outputPointSize = -pointSize / cameraSpacePos.z;
        fogDepth = pointSize / outputPointSize * 1.2;
      } else {  // Create size attenuation (if we're in 2D mode)
        const float PI = 3.1415926535897932384626433832795;
        const float minScale = 0.1;  // minimum scaling factor
        const float outSpeed = 2.0;  // shrink speed when zooming out
        const float outNorm = (1. - minScale) / atan(outSpeed);
        const float maxScale = 15.0;  // maximum scaling factor
        const float inSpeed = 0.02;  // enlarge speed when zooming in
        const float zoomOffset = 0.3;  // offset zoom pivot
        float m = projectionMatrix[0][0];
        if (zoomFactorSpecified) {
          m = zoomFactor;
        } 
        float zoom = m + zoomOffset;  // zoom pivot
        float scale = zoom < 1. ? 1. + outNorm * atan(outSpeed * (zoom - 1.)) :
                      1. + 2. / PI * (maxScale - 1.) * atan(inSpeed * (zoom - 1.));
        outputPointSize = pointSize * scale;
      }
      gl_PointSize = outputPointSize * scaleFactor;
    }`;

const FRAGMENT_SHADER_POINT_TEST_CHUNK = `
    bool point_in_unit_circle(vec2 spriteCoord) {
      vec2 centerToP = spriteCoord - vec2(0.5, 0.5);
      return dot(centerToP, centerToP) < (0.5 * 0.5);
    }

    bool point_in_unit_equilateral_triangle(vec2 spriteCoord) {
      vec3 v0 = vec3(0, 1, 0);
      vec3 v1 = vec3(0.5, 0, 0);
      vec3 v2 = vec3(1, 1, 0);
      vec3 p = vec3(spriteCoord, 0);
      float p_in_v0_v1 = cross(v1 - v0, p - v0).z;
      float p_in_v1_v2 = cross(v2 - v1, p - v1).z;
      return (p_in_v0_v1 > 0.0) && (p_in_v1_v2 > 0.0);
    }

    bool point_in_unit_square(vec2 spriteCoord) {
      return true;
    }
  `;

const FRAGMENT_SHADER = `
    
    varying vec4 vColor;

    ${ShaderChunk['common']}
    ${FRAGMENT_SHADER_POINT_TEST_CHUNK}
    uniform vec3 fogColor;
    varying float fogDepth;
		uniform float fogNear;
    uniform float fogFar;

    void main() {
      bool inside = point_in_unit_circle(gl_PointCoord);
      if (!inside) {
        discard;
      }
      gl_FragColor = vColor;
      float fogFactor = smoothstep( fogNear, fogFar, fogDepth );
      gl_FragColor.rgb = mix( gl_FragColor.rgb, fogColor, fogFactor );
    }`;


/**
 * Uses GL point sprites, either generated or from a spritesheet image to
 * render the dataset.
 */
export class ScatterPlotVisualizerSprites implements ScatterPlotVisualizer {
    public id = 'SPRITES';

    private scene!: Scene;
    private fog!: Fog;
    private renderMaterial: ShaderMaterial;
    private points!: Points;
    private worldSpacePointPositions = new Float32Array(0);
    private renderColors = new Float32Array(0);
    private standinTextureForPoints: Texture;
    private zoomFactor: number = NaN;

    constructor(private styles: Styles, spriteSheetParams?: SpriteSheetParams) {
        this.renderMaterial = this.createRenderMaterial();
        this.standinTextureForPoints = util.createTextureFromCanvas(
            document.createElement('canvas')
        );
    }

    private createUniforms(): any {
        return {
            texture: {type: 't'},
            zoomFactor: {type: 'f'},
            zoomFactorSpecified: {type: 'bool'},
            fogColor: {type: 'c'},
            fogNear: {type: 'f'},
            fogFar: {type: 'f'},
            isImage: {type: 'bool'},
            sizeAttenuation: {type: 'bool'},
            pointSize: {type: 'f'},
        };
    }

    private createRenderMaterial(): ShaderMaterial {
        const uniforms = this.createUniforms();
        return new ShaderMaterial({
            uniforms: uniforms,
            vertexShader: VERTEX_SHADER,
            fragmentShader: FRAGMENT_SHADER,
            transparent: true,
            fog: true,
            depthTest: true,
            depthWrite: true,
            blending: NormalBlending,
        });
    }


    private calculatePointSize(sceneIs3D: boolean): number {
        const n =
            this.worldSpacePointPositions != null
                ? this.worldSpacePointPositions.length / XYZ_NUM_ELEMENTS
                : 1;
        const SCALE = 200;
        const LOG_BASE = 8;
        const DIVISOR = 1.5;
        // Scale point size inverse-logarithmically to the number of points.
        const pointSize = SCALE / Math.log(n) / Math.log(LOG_BASE);
        return sceneIs3D ? pointSize : pointSize / DIVISOR;
    }

    /**
     * Create points, set their locations and actually instantiate the
     * geometry.
     */
    private createPointSprites(scene: Scene, positions: Float32Array) {
        const pointCount =
            positions != null ? positions.length / XYZ_NUM_ELEMENTS : 0;
        const geometry = this.createGeometry(pointCount);

        this.fog = new Fog(0xffffff); // unused value, gets overwritten.
        this.points = new Points(geometry, this.renderMaterial);
        this.points.frustumCulled = false;
        scene.add(this.points);
    }

    /**
     * Set up buffer attributes to be used for the points/images.
     */
    private createGeometry(pointCount: number): BufferGeometry {
        const n = pointCount;

        const geometry = new BufferGeometry();
        geometry.setAttribute(
            'position',
            new BufferAttribute(new Float32Array([]), XYZ_NUM_ELEMENTS)
        );
        geometry.setAttribute(
            'color',
            new BufferAttribute(new Float32Array([]), RGBA_NUM_ELEMENTS)
        );
        geometry.setAttribute(
            'scaleFactor',
            new BufferAttribute(new Float32Array([]), INDEX_NUM_ELEMENTS)
        );
        geometry.computeVertexNormals();
        return geometry;
    }

    private setFogDistances(
        sceneIs3D: boolean,
        nearestPointZ: number,
        farthestPointZ: number
    ) {
        const {threshold, enabled} = this.styles.fog;

        if (sceneIs3D && enabled) {
            const n = this.worldSpacePointPositions.length / XYZ_NUM_ELEMENTS;
            this.fog.near = nearestPointZ;
            // If there are fewer points we want less fog. We do this
            // by making the "far" value (that is, the distance from the camera to the
            // far edge of the fog) proportional to the number of points.
            let multiplier = 2 - Math.min(n, threshold) / threshold;
            this.fog.far = farthestPointZ * multiplier;
        } else {
            this.fog.near = Infinity;
            this.fog.far = Infinity;
        }
    }

    dispose() {
        this.disposeGeometry();
        this.disposeSpriteSheet();
    }

    private disposeGeometry() {
        if (this.points != null) {
            this.scene.remove(this.points);
            this.points.geometry.dispose();
            (this.points as any) = null;
            (this.worldSpacePointPositions as any) = null;
        }
    }

    private disposeSpriteSheet() {

        (this.renderMaterial as any) = null;
    }

    setScene(scene: Scene) {
        this.scene = scene;
    }


    onPointPositionsChanged(newPositions: Float32Array) {
        if (this.points != null) {
            if (this.worldSpacePointPositions.length !== newPositions.length) {
                this.disposeGeometry();
            }
        }

        this.worldSpacePointPositions = newPositions;
        if (this.points == null) {
            this.createPointSprites(this.scene, newPositions);
        }
        this.renderMaterial = this.createRenderMaterial();

        const positions = (this.points
            .geometry as BufferGeometry).getAttribute(
            'position'
        ) as BufferAttribute;
        positions.array = newPositions;
        positions.count = newPositions.length / XYZ_NUM_ELEMENTS;
        positions.needsUpdate = true;
    }


    onRender(rc: RenderContext) {
        const sceneIs3D: boolean = rc.camera instanceof PerspectiveCamera;
        this.setFogDistances(
            sceneIs3D,
            rc.nearestCameraSpacePointZ,
            rc.farthestCameraSpacePointZ
        );
        this.scene.fog = this.fog;
        this.scene.fog.color = new Color(rc.backgroundColor);
        this.renderMaterial.uniforms.fogColor.value = this.scene.fog.color;
        this.renderMaterial.uniforms.fogNear.value = this.fog.near;
        this.renderMaterial.uniforms.fogFar.value = this.fog.far;
        this.renderMaterial.uniforms.texture.value = this.standinTextureForPoints;
        this.renderMaterial.uniforms.sizeAttenuation.value = sceneIs3D;
        this.renderMaterial.uniforms.zoomFactor.value = this.zoomFactor;
        this.renderMaterial.uniforms.zoomFactorSpecified.value = !isNaN(this.zoomFactor);
        this.renderMaterial.uniforms.pointSize.value = this.calculatePointSize(
            sceneIs3D
        );
        this.points.material = this.renderMaterial;
        let colors = (this.points.geometry as BufferGeometry).getAttribute(
            'color'
        ) as BufferAttribute;
        this.renderColors = rc.pointColors;
        colors.array = this.renderColors;
        colors.count = this.renderColors.length / RGBA_NUM_ELEMENTS;
        colors.needsUpdate = true;
        let scaleFactors = (this.points
            .geometry as BufferGeometry).getAttribute(
            'scaleFactor'
        ) as BufferAttribute;
        scaleFactors.array = rc.pointScaleFactors;
        scaleFactors.count = rc.pointScaleFactors.length / INDEX_NUM_ELEMENTS;
        scaleFactors.needsUpdate = true;
    }

    onResize(newWidth: number, newHeight: number) {
    }
}

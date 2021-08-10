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
import {Styles} from './styles';
import * as util from './util';
import {RGB_NUM_ELEMENTS, RGBA_NUM_ELEMENTS, UV_NUM_ELEMENTS, XYZ_NUM_ELEMENTS,} from './constants';
import {BufferAttribute, BufferGeometry, Color, Mesh, Scene, ShaderMaterial, Texture} from "three";

const MAX_CANVAS_DIMENSION = 8192;
const NUM_GLYPHS = 256;
const VERTICES_PER_GLYPH = 2 * 3; // 2 triangles, 3 verts per triangle

/**
 * Each label is made up of triangles (two per letter.) Each vertex, then, is
 * the corner of one of these triangles (and thus the corner of a letter
 * rectangle.)
 * Each has the following attributes:
 *    posObj: The (x, y) position of the vertex within the label, where the
 *            bottom center of the word is positioned at (0, 0);
 *    position: The position of the label in worldspace.
 *    vUv: The (u, v) coordinates that index into the glyphs sheet (range 0, 1.)
 *    color: The color of the label (matches the corresponding point's color.)
 *    wordShown: Boolean. Whether or not the label is visible.
 */

const makeVertexShader = (fontSize: number, scale: number) => `
      attribute vec2 posObj;
      attribute vec4 color;
      varying vec2 vUv;
      varying vec4 vColor;

      void main() {
        vUv = uv;
        vColor = color;

        // Rotate label to face camera.

        vec4 vRight = vec4(
          modelViewMatrix[0][0], modelViewMatrix[1][0], modelViewMatrix[2][0], 0);

        vec4 vUp = vec4(
          modelViewMatrix[0][1], modelViewMatrix[1][1], modelViewMatrix[2][1], 0);

        vec4 vAt = -vec4(
          modelViewMatrix[0][2], modelViewMatrix[1][2], modelViewMatrix[2][2], 0);

        mat4 pointToCamera = mat4(vRight, vUp, vAt, vec4(0, 0, 0, 1));

        vec2 scaledPos = posObj * ${1 / fontSize} * ${scale};

        vec4 posRotated = pointToCamera * vec4(scaledPos, 0, 1);
        vec4 mvPosition = modelViewMatrix * (vec4(position, 0) + posRotated);
        gl_Position = projectionMatrix * mvPosition;
      }`;

const FRAGMENT_SHADER = `
      uniform sampler2D texture;
      // uniform bool picking;
      varying vec2 vUv;
      varying vec4 vColor;

      void main() {
        // if (picking) {
        //   gl_FragColor = vColor;
        // } else {
          vec4 fromTexture = texture2D(texture, vUv);
          gl_FragColor = vColor * fromTexture;
        // }
      }`;

type GlyphTexture = {
  texture: Texture;
  lengths: Float32Array;
  offsets: Float32Array;
};

/**
 * Renders the text labels as 3d geometry in the world.
 */
export class ScatterPlotVisualizer3DLabels implements ScatterPlotVisualizer {
  public id = '3D_LABELS';

  private scene!: Scene;
  private labelStrings: string[] = [];
  private geometry!: BufferGeometry;
  private worldSpacePointPositions = new Float32Array(0);
  // private pickingColors = new Float32Array(0);
  private renderColors = new Float32Array(0);
  private material!: ShaderMaterial;
  private uniforms: any = {};
  private labelsMesh!: Mesh;
  private positions!: BufferAttribute;
  private totalVertexCount = 0;
  private labelVertexMap: number[][] = [];
  private glyphTexture!: GlyphTexture;

  constructor(private styles: Styles) {}

  private createGlyphTexture(): GlyphTexture {
    const {fontSize, backgroundColor, color} = this.styles.label3D;
    const canvas = document.createElement('canvas');
    canvas.width = MAX_CANVAS_DIMENSION;
    canvas.height = fontSize;
    const ctx = canvas.getContext('2d')!;
    ctx.font = 'bold ' + fontSize + 'px Roboto Condensed';
    ctx.textBaseline = 'top';
    ctx.fillStyle = backgroundColor;
    // ctx.rect(0, 0, canvas.width, canvas.height);
    // ctx.fill();

    ctx.fillStyle = color;
    const spaceOffset = ctx.measureText(' ').width;
    // For each letter, store length, position at the encoded index.
    const glyphLengths = new Float32Array(NUM_GLYPHS);
    const glyphOffset = new Float32Array(NUM_GLYPHS);
    let leftCoord = 0;
    for (let i = 0; i < NUM_GLYPHS; i++) {
      const text = ' ' + String.fromCharCode(i);
      const textLength = ctx.measureText(text).width;
      glyphLengths[i] = textLength - spaceOffset;
      glyphOffset[i] = leftCoord;
      ctx.fillText(text, leftCoord - spaceOffset, 0);
      leftCoord += textLength;
    }
    const tex = util.createTextureFromCanvas(canvas);
    return {texture: tex, lengths: glyphLengths, offsets: glyphOffset};
  }

  private processLabelVerts(pointCount: number) {
    let numTotalLetters = 0;
    this.labelVertexMap = [];
    for (let i = 0; i < pointCount; i++) {
      const label = this.labelStrings[i];
      let vertsArray: number[] = [];
      for (let j = 0; j < label.length; j++) {
        for (let k = 0; k < VERTICES_PER_GLYPH; k++) {
          vertsArray.push(numTotalLetters * VERTICES_PER_GLYPH + k);
        }
        numTotalLetters++;
      }
      this.labelVertexMap.push(vertsArray);
    }
    this.totalVertexCount = numTotalLetters * VERTICES_PER_GLYPH;
  }

  private createColorBuffers(pointCount: number) {
    // this.pickingColors = new Float32Array(
    //   this.totalVertexCount * RGB_NUM_ELEMENTS
    // );
    this.renderColors = new Float32Array(
      this.totalVertexCount * RGB_NUM_ELEMENTS
    );
    for (let i = 0; i < pointCount; i++) {
      // const pickingColor = new Color(i);
      this.labelVertexMap[i].forEach(j => {
        // this.pickingColors[RGB_NUM_ELEMENTS * j] = pickingColor.r;
        // this.pickingColors[RGB_NUM_ELEMENTS * j + 1] = pickingColor.g;
        // this.pickingColors[RGB_NUM_ELEMENTS * j + 2] = pickingColor.b;
        this.renderColors[RGB_NUM_ELEMENTS * j] = 1.0;
        this.renderColors[RGB_NUM_ELEMENTS * j + 1] = 1.0;
        this.renderColors[RGB_NUM_ELEMENTS * j + 2] = 1.0;

      });
    }
  }

  private createLabels() {
    const {fontSize, scale} = this.styles.label3D;
    if (this.labelStrings == null || this.worldSpacePointPositions == null) {
      return;
    }
    const pointCount = this.worldSpacePointPositions.length / XYZ_NUM_ELEMENTS;
    if (pointCount !== this.labelStrings.length) {
      console.log(pointCount + '!=' +  this.labelStrings.length)
      return;
    }
    this.glyphTexture = this.createGlyphTexture();

    this.uniforms = {
      texture: {type: 't'},
      picking: {type: 'bool'},
    };


    this.material = new ShaderMaterial({
      uniforms: this.uniforms,
      transparent: true,
      vertexShader: makeVertexShader(fontSize, scale),
      fragmentShader: FRAGMENT_SHADER,
      // polygonOffset: true,
      // polygonOffsetFactor: -4,
      depthWrite: false,
      depthTest:false
    });

    this.processLabelVerts(pointCount);
    this.createColorBuffers(pointCount);

    let positionArray = new Float32Array(
      this.totalVertexCount * XYZ_NUM_ELEMENTS
    );
    this.positions = new BufferAttribute(positionArray, XYZ_NUM_ELEMENTS);

    let posArray = new Float32Array(this.totalVertexCount * XYZ_NUM_ELEMENTS);
    let uvArray = new Float32Array(this.totalVertexCount * UV_NUM_ELEMENTS);
    let colorsArray = new Float32Array(
      this.totalVertexCount * RGB_NUM_ELEMENTS
    );
    let positionObject = new BufferAttribute(posArray, 2);
    let uv = new BufferAttribute(uvArray, UV_NUM_ELEMENTS);
    let colors = new BufferAttribute(colorsArray, RGB_NUM_ELEMENTS);

    this.geometry = new BufferGeometry();
    this.geometry.setAttribute('posObj', positionObject);
    this.geometry.setAttribute('position', this.positions);
    this.geometry.setAttribute('uv', uv);
    this.geometry.setAttribute('color', colors);

    let lettersSoFar = 0;
    for (let i = 0; i < pointCount; i++) {
      const label = this.labelStrings[i];
      let leftOffset = 0;
      // Determine length of word in pixels.
      for (let j = 0; j < label.length; j++) {
        let letterCode = label.charCodeAt(j);
        leftOffset += this.glyphTexture.lengths[letterCode];
      }
      leftOffset /= -2; // centers text horizontally around the origin
      for (let j = 0; j < label.length; j++) {
        let letterCode = label.charCodeAt(j);
        let letterWidth = this.glyphTexture.lengths[letterCode];
        let scale = fontSize;
        let right = (leftOffset + letterWidth) / scale;
        let left = leftOffset / scale;
        let top = fontSize / scale;

        // First triangle
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 0, left, 0);
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 1, right, 0);
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 2, left, top);

        // Second triangle
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 3, left, top);
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 4, right, 0);
        positionObject.setXY(lettersSoFar * VERTICES_PER_GLYPH + 5, right, top);

        // Set UVs based on letter.
        let uLeft = this.glyphTexture.offsets[letterCode];
        let uRight = this.glyphTexture.offsets[letterCode] + letterWidth;
        // Scale so that uvs lie between 0 and 1 on the texture.
        uLeft /= MAX_CANVAS_DIMENSION;
        uRight /= MAX_CANVAS_DIMENSION;
        let vTop = 1;
        let vBottom = 0;
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 0, uLeft, vTop);
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 1, uRight, vTop);
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 2, uLeft, vBottom);
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 3, uLeft, vBottom);
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 4, uRight, vTop);
        uv.setXY(lettersSoFar * VERTICES_PER_GLYPH + 5, uRight, vBottom);

        lettersSoFar++;
        leftOffset += letterWidth;
      }
    }

    for (let i = 0; i < pointCount; i++) {
      const p = util.vector3FromPackedArray(this.worldSpacePointPositions, i);
      this.labelVertexMap[i].forEach(j => {
        this.positions.setXYZ(j, p.x, p.y, p.z);
      });
    }

    this.labelsMesh = new Mesh(this.geometry, this.material);
    this.labelsMesh.name = '3d_labels_text';
    this.labelsMesh.renderOrder = 999;
    this.labelsMesh.frustumCulled = false;

  }


  private colorLabels() {
    if (
      this.labelStrings == null ||
      this.geometry == null) {
      return;
    }
    const colors = this.geometry.getAttribute('color') as BufferAttribute;
    colors.array = this.renderColors;

    const n = this.labelStrings.length;
    let src = 0;
    const c = new Color(
        1,
        1,
        1
    );
    for (let i = 0; i < n; ++i) {
      // const c = new Color(
      //   pointColors[src],
      //   pointColors[src + 1],
      //   pointColors[src + 2]
      // );

      const m = this.labelVertexMap[i].length;
      for (let j = 0; j < m; ++j) {
        colors.setXYZ(this.labelVertexMap[i][j], c.r, c.g, c.b);
      }
      src += RGBA_NUM_ELEMENTS;
    }
    colors.needsUpdate = true;
  }

  setScene(scene: Scene) {
    this.scene = scene;
  }

  dispose() {
    if (this.labelsMesh) {
      if (this.scene) {
        this.scene.remove(this.labelsMesh);
      }
      (this.labelsMesh as any) = null;
    }
    if (this.geometry) {
      this.geometry.dispose();
      (this.geometry as any) = null;
    }
    if (this.glyphTexture != null && this.glyphTexture.texture != null) {
      this.glyphTexture.texture.dispose();
      (this.glyphTexture as any).texture = null;
    }
  }

  onPickingRender(rc: RenderContext) {
    // if (this.geometry == null) {
    //   this.createLabels();
    //   return;
    // }
    // this.material.uniforms.texture.value = this.glyphTexture.texture;
    // this.material.uniforms.picking.value = true;
    // const colors = this.geometry.getAttribute('color') as BufferAttribute;
    // colors.array = this.pickingColors;
    // colors.needsUpdate = true;
    this.removeMesh();
  }

  private removeMesh() {
    this.scene.remove(this.labelsMesh);
  }
  onRender(rc: RenderContext) {
    if (this.geometry == null) {
      this.createLabels();

    }
    this.removeMesh();
    this.scene.add(this.labelsMesh);
    this.colorLabels();
    this.material.uniforms.texture.value = this.glyphTexture.texture;
    this.material.uniforms.picking.value = false;
    const colors = this.geometry.getAttribute('color') as BufferAttribute;
    colors.array = this.renderColors;
    colors.needsUpdate = true;
  }

  onPointPositionsChanged(newPositions: Float32Array) {
    // ignore
  }

  setLabelStrings(labelStrings: string[]) {
    this.labelStrings = labelStrings;
    this.dispose();
  }

  setLabels(labelStrings: string[], positions: Float32Array) {
    this.labelStrings = labelStrings;
    this.worldSpacePointPositions = positions;
    this.dispose();
  }

  onResize(newWidth: number, newHeight: number) {}
}

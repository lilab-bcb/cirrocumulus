import Link from '@material-ui/core/Link';
import React from 'react';

function stateReset(gl) {
    let numAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
    let tmp = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, tmp);
    for (let ii = 0; ii < numAttribs; ++ii) {
        gl.disableVertexAttribArray(ii);
        gl.vertexAttribPointer(ii, 4, gl.FLOAT, false, 0, 0);
        gl.vertexAttrib1f(ii, 0);
    }
    gl.deleteBuffer(tmp);

    let numTextureUnits = gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS);
    for (let ii = 0; ii < numTextureUnits; ++ii) {
        gl.activeTexture(gl.TEXTURE0 + ii);
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, null);
        gl.bindTexture(gl.TEXTURE_2D, null);
    }

    gl.activeTexture(gl.TEXTURE0);
    gl.useProgram(null);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.bindRenderbuffer(gl.RENDERBUFFER, null);
    gl.disable(gl.BLEND);
    gl.disable(gl.CULL_FACE);
    gl.disable(gl.DEPTH_TEST);
    gl.disable(gl.DITHER);
    gl.disable(gl.SCISSOR_TEST);
    gl.blendColor(0, 0, 0, 0);
    gl.blendEquation(gl.FUNC_ADD);
    gl.blendFunc(gl.ONE, gl.ZERO);
    gl.clearColor(0, 0, 0, 0);
    gl.clearDepth(1);
    gl.clearStencil(-1);
    gl.colorMask(true, true, true, true);
    gl.cullFace(gl.BACK);
    gl.depthFunc(gl.LESS);
    gl.depthMask(true);
    gl.depthRange(0, 1);
    gl.frontFace(gl.CCW);
    gl.hint(gl.GENERATE_MIPMAP_HINT, gl.DONT_CARE);
    gl.lineWidth(1);
    gl.pixelStorei(gl.PACK_ALIGNMENT, 4);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
    gl.polygonOffset(0, 0);
    gl.sampleCoverage(1, false);
    gl.scissor(0, 0, gl.canvas.width, gl.canvas.height);
    gl.stencilFunc(gl.ALWAYS, 0, 0xFFFFFFFF);
    gl.stencilMask(0xFFFFFFFF);
    gl.stencilOp(gl.KEEP, gl.KEEP, gl.KEEP);
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT | gl.STENCIL_BUFFER_BIT);

    return gl;
}

class GalleryImage extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {url: null};
    }


    componentDidUpdate(prevProps, prevState) {

        if (this.state.url == null) {

            const {traceInfo, config} = this.props;
            const div = document.createElement('div');
            const layout = Object.assign({}, traceInfo.layout);
            const myTraceInfo = Object.assign({}, traceInfo);
            myTraceInfo.data.forEach(trace => {
                trace.selectedpoints = null;
                // do not draw selection
                // TODO 3d plots
            });

            layout.width = 400;
            layout.height = 400;

            window.Plotly.react(div, {
                data: traceInfo.data,
                layout: layout,
                config: config
            });
            // "gl-canvas gl-canvas-context"
            // "gl-canvas gl-canvas-focus"
            // "gl-canvas gl-canvas-pick"
            let canvas = div.querySelectorAll('.gl-canvas');
            window.Plotly.purge(div);
            this.setState({url: canvas[0].toDataURL()});

            // div.querySelectorAll('canvas').forEach(canvas => {
            //     console.log('here');
            //     stateReset(canvas.getContext('webgl'));
            // });

        }
    }


    onSelect = (event) => {
        event.preventDefault();
        this.props.onSelect(this.props.traceInfo);
    };

    render() {
        if (this.state.url == null) {
            return null;
        }
        return (
            <div style={{display: 'inline-block'}}>
                <Link href="#"
                      onClick={this.onSelect}>{this.props.traceInfo.data[0].name}</Link>
                {/*<div style={{position: 'relative'}}>*/}
                {/*    {this.state.urls.map(url =>*/}
                {/*        <img style={{position: 'absolute', top: 0, left: 0}} src={url}/>)}*/}
                {/*</div>*/}
                <img src={this.state.url}/>
            </div>
        );
    }
}

export default GalleryImage;


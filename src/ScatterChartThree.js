import React from 'react';
import {Dataset, ScatterGL} from 'scatter-gl';


import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {drawScatter2d, getChartSize} from './PlotUtil';

class ScatterChartThree extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();
        this.scatterGL = null;
        this.chartSize = getChartSize();
        this.state = {animating: false};

        // window.addEventListener('resize', () => {
        //     scatterGL.resize();
        // });

    }


    static snapshot(scatterGL, traceInfo, markerOpacity) {
        const dataset = new Dataset(traceInfo.x, traceInfo.y, traceInfo.z, traceInfo.marker.color);
        scatterGL.setPointColorer((i, selectedIndices, hoverIndex) => {
            const c = dataset.metadata[i];
            c.opacity = markerOpacity;
            return c;
        });
        scatterGL.setDimensions(traceInfo.z == null ? 2 : 3);
        scatterGL.render(dataset);
        // scatterGL.updateScatterPlotAttributes();
        // scatterGL.renderScatterPlot();
    }


    componentDidMount() {
        this.init();
        this.draw();
    }

    componentWillUnmount() {
        console.log('unmount');
        // if (this.scatterGL != null) {
        //     this.scatterGL.dispose();
        // }
    }


    onHome = () => {

    };

    onSaveImage = () => {
        const markerSize = 4; //this.scatterGL.pointVisualizer.renderMaterial.uniforms.pointSize.value;
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, color} = this.props;
        const chartSize = {width: 800, height: 800};
        let context = new window.C2S(chartSize.width, chartSize.height);
        drawScatter2d(context, chartSize, traceInfo, markerSize, markerOpacity, unselectedMarkerOpacity, selection, color);
        let svg = context.getSerializedSvg();
        let blob = new Blob([svg], {
            type: 'text/plain;charset=utf-8'
        });
        let name = traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        // const renderer = new SVGRenderer();
        // document.body.appendChild(renderer.domElement);
        // renderer.setSize(window.innerWidth, window.innerHeight);
        // renderer.setQuality('low');
        // renderer.render(this.scatterGL.scatterPlot.scene, this.scatterGL.scatterPlot.camera);
        // console.log(renderer.domElement);
        window.saveAs(blob, name + '.svg');
    };
    onToggleAnimation = () => {
        if (this.scatterGL.scatterPlot.orbitIsAnimating()) {
            this.scatterGL.stopOrbitAnimation();
            this.setState({animating: false});
        } else {
            this.scatterGL.startOrbitAnimation();
            this.setState({animating: true});
        }
    };

    onDragMode = (mode) => {
        if (mode === 'pan') {
            this.scatterGL.setPanMode();
        } else if (mode === 'select') {
            this.scatterGL.setSelectMode();
        }
    };

    init() {
        const {traceInfo} = this.props;
        if (this.scatterGL == null) {
            this.scatterGL = new ScatterGL(this.containerElementRef.current, {
                renderMode: 'POINT',
                rotateOnStart: false,
                showLabelsOnHover: false,
                onSelect: (selectedpoints, boundingBox) => {
                    if (this.scatterGL.scatterPlot.interactionMode === 'PAN') {
                        return;
                    }
                    if (selectedpoints != null && selectedpoints.length === 0) {
                        selectedpoints = null;
                    }

                    if (selectedpoints == null) {
                        this.props.onDeselect({name: getEmbeddingKey(traceInfo.embedding)});
                    } else {

                        let xmin = Number.MAX_VALUE;
                        let ymin = Number.MAX_VALUE;
                        let zmin = Number.MAX_VALUE;
                        let xmax = -Number.MAX_VALUE;
                        let ymax = -Number.MAX_VALUE;
                        let zmax = -Number.MAX_VALUE;
                        const is3d = traceInfo.z != null;
                        selectedpoints.forEach(index => {
                            const x = traceInfo.x[index];
                            xmin = Math.min(xmin, x);
                            xmax = Math.max(xmax, x);
                            const y = traceInfo.y[index];
                            ymin = Math.min(ymin, y);
                            ymax = Math.max(ymax, y);
                            if (is3d) {
                                const z = traceInfo.z[index];
                                zmin = Math.min(zmin, z);
                                zmax = Math.max(zmax, z);
                            }
                        });


                        let path = {shape: 'rect', x: xmin, y: ymin, width: xmax - xmin, height: ymax - ymin};
                        if (is3d) {
                            path.shape = 'rect 3d';
                            path.z = zmin;
                            path.depth = zmax - zmin;
                        }
                        this.props.onSelected({
                            name: getEmbeddingKey(traceInfo.embedding),
                            value: {basis: traceInfo.embedding, selectedpoints: selectedpoints, path: path}
                        });
                    }
                },
                onHover: (point) => {


                }
            });

            this.scatterGL.setSelectMode();
            const canvas = this.containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
        }
    }

    draw() {
        const scatterGL = this.scatterGL;
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, color} = this.props;
        scatterGL.setSelectedPointIndices(selection);
        const dataset = new Dataset(traceInfo.x, traceInfo.y, traceInfo.z, color);
        scatterGL.render(dataset);
        scatterGL.setDimensions(traceInfo.z == null ? 2 : 3);
        scatterGL.setPointColorer((i, selectedIndices, hoverIndex) => {
            const c = dataset.metadata[i];
            c.opacity = markerOpacity;
            // if (hoverIndex === i) {
            //     return c.brighter();
            // }
            const isSelected = selectedIndices.size === 0 || selectedIndices.has(i);
            if (!isSelected) {
                c.opacity = unselectedMarkerOpacity;
            }
            return c;
        });
        scatterGL.render(dataset);
        scatterGL.updateScatterPlotAttributes();
        scatterGL.renderScatterPlot();
        window.foo = this.scatterGL;
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        this.draw();
    }

    render() {

        return <React.Fragment><ChartToolbar onHome={this.onHome}
                                             is3d={this.props.traceInfo && this.props.traceInfo.z != null}
                                             animating={this.state.animating}
                                             toggleAnimation={this.onToggleAnimation}
                                             onSaveImage={this.onSaveImage}
                                             onDragMode={this.onDragMode}></ChartToolbar>
            <div style={{
                display: 'inline-block',
                position: 'relative',
                width: this.chartSize.width,
                height: this.chartSize.height
            }}
                 ref={this.containerElementRef}/>
        </React.Fragment>;
    }
}

export default ScatterChartThree;




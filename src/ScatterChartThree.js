import Link from '@material-ui/core/Link';
import React from 'react';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {drawScatter2d, getChartSize} from './PlotUtil';
import {createScatterPlot} from './ThreeUtil';

class ScatterChartThree extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();
        this.scatterPlot = null;
        this.chartSize = getChartSize();
        this.state = {animating: false};

        // window.addEventListener('resize', () => {
        //     scatterGL.resize();
        // });

    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.draw();

    }

    componentDidMount() {
        this.init();
        this.draw();
        this.draw();
    }

    componentWillUnmount() {
        console.log('unmount');
        // if (this.scatterGL != null) {
        //     this.scatterGL.dispose();
        // }
    }


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
        if (this.scatterPlot.orbitIsAnimating()) {
            this.scatterPlot.stopOrbitAnimation();
            this.setState({animating: false});
        } else {
            this.scatterPlot.startOrbitAnimation();
            this.setState({animating: true});
        }
    };

    onDragMode = (mode) => {
        if (mode === 'pan') {
            this.scatterPlot.setInteractionMode('PAN');
        } else if (mode === 'select') {
            this.scatterPlot.setInteractionMode('SELECT');
        }
    };

    onGallery = (event) => {
        event.preventDefault();
        this.props.onGallery();
    };

    init() {
        const {traceInfo} = this.props;
        if (this.scatterPlot == null) {
            const containerElement = this.containerElementRef.current;
            this.scatterPlot = createScatterPlot(containerElement);
            this.scatterPlot.selectCallback = (selectedpoints) => {
                if (this.scatterPlot.interactionMode === 'PAN') {
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
            };

            // this.scatterGL = new ScatterGL(this.containerElementRef.current, {
            //     renderMode: 'POINT',
            //     rotateOnStart: false,
            //     showLabelsOnHover: false,
            //
            //     onHover: (point) => {
            //
            //     }
            // });


            const canvas = this.containerElementRef.current.querySelector('canvas');
            canvas.style.outline = '0px';
        }
    }

    draw() {
        const {traceInfo, markerOpacity, unselectedMarkerOpacity, selection, color, pointSize} = this.props;
        const colors = traceInfo.colors;
        for (let i = 0, j = 3, n = traceInfo.npoints; i < n; i++, j += 4) {
            const isSelected = selection.size === 0 || selection.has(i);
            colors[j] = isSelected ? markerOpacity : unselectedMarkerOpacity;
        }
        this.scatterPlot.setPointColors(colors);
        this.scatterPlot.setPointPositions(traceInfo.positions);
        this.scatterPlot.setDimensions(traceInfo.dimensions);
        // const {scaleDefault, scaleSelected, scaleHover} = this.scatterPlot.styles.point;

        const scale = new Float32Array(traceInfo.npoints);
        scale.fill(pointSize);
        this.scatterPlot.setPointScaleFactors(scale);
        this.scatterPlot.render();
    }


    render() {
        return <React.Fragment><ChartToolbar onHome={this.onHome}
                                             is3d={this.props.traceInfo && this.props.traceInfo.z != null}
                                             animating={this.state.animating}
                                             toggleAnimation={this.onToggleAnimation}
                                             onSaveImage={this.onSaveImage}
                                             onDragMode={this.onDragMode}></ChartToolbar>
            <Link style={{paddingLeft: 5}} href="#" onClick={this.onGallery}>
                Gallery
            </Link>
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




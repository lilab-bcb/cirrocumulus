import withStyles from '@material-ui/core/styles/withStyles';
import Typography from '@material-ui/core/Typography';
import React from 'react';
import svgPanZoom from 'svg-pan-zoom';
import ChartToolbar from './ChartToolbar';
import {intFormat, numberFormat2f} from './formatters';

const styles = theme => ({

    root: {
        '& > *': {
            margin: theme.spacing(.4),
        },
        '& > .MuiIconButton-root': {
            padding: 0,
        },
        position: 'absolute',
        zIndex: 1,
        top: 0,
        left: 0,
        display: 'inline-block',
        verticalAlign: 'top',
        whiteSpace: 'nowrap',
        overflow: 'hidden',
    }
});

class MetaEmbedding extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();

        this.tooltipElementRef = React.createRef();
        this.state = {loading: false};
    }


    onZoomIn = () => {
        this.zoom.zoomIn();
    };

    onZoomOut = () => {
        this.zoom.zoomOut();
    };


    onSaveImage = (format) => {
        let context;
        let canvas = null;
        const {chartSize, traceInfo} = this.props;
        const totalSize = {width: chartSize.width, height: chartSize.height};
        let name = traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        if (format !== 'svg') {
            canvas = document.createElement('canvas');
            canvas.width = totalSize.width * window.devicePixelRatio;
            canvas.height = totalSize.height * window.devicePixelRatio;
            context = canvas.getContext('2d');
            context.scale(window.devicePixelRatio, window.devicePixelRatio);
            context.fillStyle = 'white';
            context.fillRect(0, 0, totalSize.width, totalSize.height);
            const xml = new XMLSerializer().serializeToString(traceInfo.source);
            const svg64 = btoa(xml);
            const b64Start = 'data:image/svg+xml;base64,';
            const image64 = b64Start + svg64;
            context.fillStyle = 'black';
            const img = new Image();
            img.src = image64;
            img.onload = function () {
                context.drawImage(img, 0, 0);
                canvas.toBlob(blob => {
                    window.saveAs(blob, name + '.png', true);
                });
            };
        } else {
            let blob = new Blob([new XMLSerializer().serializeToString(traceInfo.source)], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, name + '.svg');
        }
    };

    onHome = () => {
        this.zoom.fit();
        this.zoom.center();
    };
    updateSvg = () => {
        const containerElement = this.containerElementRef.current;
        containerElement.innerHTML = '';
        const svg = this.props.traceInfo.source;
        svg.setAttribute('width', this.props.chartSize.width);
        svg.setAttribute('height', this.props.chartSize.height);
        containerElement.append(svg);
        this.zoom = svgPanZoom(svg, {dblClickZoomEnabled: false, contain: true});
        this.zoom.fit();
    };

    componentDidMount() {
        this.updateSvg();
        const containerElement = this.containerElementRef.current;
        containerElement.addEventListener('mousemove', (e) => {
            if (e.target.nodeName === 'path') {
                const categoryToStats = this.props.traceInfo.categoryToStats;
                let text = e.target.id;
                text = text.replaceAll('_', ' '); // FIXME
                const stats = categoryToStats[text];
                if (stats) {
                    if (this.props.traceInfo.name !== '__count') {
                        if (this.props.traceInfo.isCategorical) {
                            text += ', mode: ' + stats.value + ', # spots: ' + intFormat(stats.n);
                        } else {
                            text += ', mean: ' + numberFormat2f(stats.value) + ', # spots: ' + intFormat(stats.n);
                        }
                    } else {
                        text += ', # spots: ' + intFormat(stats.n);
                    }
                }
                this.tooltipElementRef.current.innerHTML = text;
            } else {
                this.tooltipElementRef.current.innerHTML = '';
            }
        });
        containerElement.addEventListener('mouseleave', (e) => {
            this.tooltipElementRef.current.innerHTML = '';
        });
    }

    componentDidUpdate(prevProps, prevState, snapshot) {
        if (this.props.traceInfo.source !== prevProps.traceInfo.source) {
            this.updateSvg();
        } else if (this.props.chartSize.width !== prevProps.chartSize.width || this.props.chartSize.height !== prevProps.chartSize.height) {
            const svg = this.props.traceInfo.source;
            svg.setAttribute('width', this.props.chartSize.width);
            svg.setAttribute('height', this.props.chartSize.height);
            this.zoom.fit();
        }
    }

    render() {

        return <React.Fragment>
            <div className={this.props.classes.root}>
                <ChartToolbar
                    // dragmode={this.props.chartOptions.dragmode}
                    // editSelection={this.props.chartOptions.editSelection}
                    // onMoreOptions={this.props.onMoreOptions}
                    onGallery={this.props.onGallery}
                    animating={false}
                    onZoomIn={this.onZoomIn}
                    onZoomOut={this.onZoomOut}
                    is3d={false}
                    onHome={this.onHome}
                    onSaveImage={this.onSaveImage}
                    // onDragMode={this.onDragMode}
                    // onEditSelection={this.onEditSelection}
                >
                </ChartToolbar>
                <Typography color="textPrimary" ref={this.tooltipElementRef} style={{
                    display: 'inline-block',
                    paddingLeft: 5,
                    verticalAlign: 'top'
                }}>&nbsp;</Typography>
            </div>

            <div style={{
                display: 'inline-block',
                // width: this.props.chartSize.width,
                // height: this.props.chartSize.height
            }}
                 ref={this.containerElementRef}>
            </div>
        </React.Fragment>;
    }
}


export default withStyles(styles)(MetaEmbedding);

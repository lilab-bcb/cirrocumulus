import withStyles from '@mui/styles/withStyles';
import {select} from 'd3-selection';
import {zoom, zoomIdentity} from 'd3-zoom';
import React from 'react';
import {getEmbeddingKey} from './actions';
import ChartToolbar from './ChartToolbar';
import {intFormat, numberFormat2f} from './formatters';
import {stripTrailingZeros} from './util';

const styles = theme => ({

    root: {
        '& > *': {
            margin: theme.spacing(.4)
        },
        '& > .MuiIconButton-root': {
            padding: 0
        },
        position: 'absolute',
        zIndex: 1,
        top: 0,
        left: 0,
        display: 'inline-block',
        verticalAlign: 'top',
        whiteSpace: 'nowrap',
        overflow: 'hidden'
    }
});

export function createCategoryToStats(trace, selection) {
    const categoryToStats = {};
    const categoryToIndices = trace.categoryToIndices;
    const selectionEmpty = selection == null;
    for (const category in categoryToIndices) {
        const indices = categoryToIndices[category];
        const valueToCount = {};
        let sum = 0;
        let n = 0;
        for (let i = 0, nIndices = indices.length; i < nIndices; i++) {
            const index = indices[i];
            if (selectionEmpty || selection.has(index)) {
                const val = trace.values[index];
                if (!trace.continuous) {
                    valueToCount[val] = (valueToCount[val] || 0) + 1;
                } else {
                    sum += val;
                }
                n++;
            }
        }
        if (!trace.continuous) {
            let maxCount = 0;
            let maxValue;
            for (let value in valueToCount) {
                const count = valueToCount[value];
                if (count > maxCount) {
                    maxCount = count;
                    maxValue = value;
                }
            }
            categoryToStats[category] = {value: maxValue, n: n};
        } else {
            const mean = sum / n;
            categoryToStats[category] = {value: n === 0 ? undefined : (mean - trace.mean) / trace.stdev, n: n};
        }
    }
    return categoryToStats;
}

class MetaEmbedding extends React.PureComponent {

    constructor(props) {
        super(props);
        this.containerElementRef = React.createRef();
    }

    onSaveImage = (format) => {
        let context;
        let canvas = null;
        const {chartSize, trace} = this.props;
        const totalSize = {width: chartSize.width, height: chartSize.height};
        let name = trace.name;
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
            const xml = new XMLSerializer().serializeToString(trace.source);
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
            let blob = new Blob([new XMLSerializer().serializeToString(trace.source)], {
                type: 'text/plain;charset=utf-8'
            });
            window.saveAs(blob, name + '.svg');
        }
    };

    updateSvg = () => {
        const containerElement = this.containerElementRef.current;
        containerElement.innerHTML = '';
        const svg = this.props.trace.source;
        svg.setAttribute('width', this.props.chartSize.width);
        svg.setAttribute('height', this.props.chartSize.height);
        containerElement.append(svg);
        let g;
        let childNodes = svg.childNodes || svg.children;
        if (childNodes.length > 0 && childNodes[0].nodeName === "g" && childNodes[0].getAttribute("cirro-zoom") === 'true') {
            g = childNodes[0];
        }
        if (!g) {
            const nodeNamesKeep = new Set(['defs', 'metadata', 'style']);
            g = document.createElementNS("http://www.w3.org/2000/svg", "g");
            g.setAttribute('cirro-zoom', 'true');
            const svgChildren = svg.childNodes || svg.children;
            if (!!svgChildren && svgChildren.length > 0) {
                for (let i = svgChildren.length; i > 0; i--) {
                    if (!nodeNamesKeep.has(svgChildren[svgChildren.length - i].nodeName)) {
                        g.appendChild(svgChildren[svgChildren.length - i]);
                    }
                }
            }
            svg.appendChild(g);
        }

        function zoomed({transform}) {
            g.setAttribute("transform", transform);
        }


        const d3Zoom = zoom().scaleExtent([0.1, 10]).on("zoom", zoomed);
        this.d3Zoom = d3Zoom;
        select(svg).call(d3Zoom);

    };

    onHome = () => {
        select(this.props.trace.source).call(this.d3Zoom.transform, zoomIdentity.scale(1));
    };

    onZoomIn = () => {
        select(this.props.trace.source).call(this.d3Zoom.scaleBy, 1.5);
    };

    onZoomOut = () => {
        select(this.props.trace.source).call(this.d3Zoom.scaleBy, .5);
    };

    componentDidMount() {
        this.updateSvg();
        const containerElement = this.containerElementRef.current;
        containerElement.addEventListener('mousemove', (e) => {
            containerElement.style.cursor = null;
            if (e.target.nodeName === 'path') {
                const categoryToStats = this.props.trace.categoryToStats;
                let category = e.target.id;
                category = category.replaceAll('_', ' '); // FIXME
                let tooltip = category;
                if (categoryToStats) {
                    containerElement.style.cursor = 'pointer';
                    const stats = categoryToStats[category];
                    if (stats) {
                        const fullStats = this.props.trace.fullCategoryToStats[category];
                        const showFull = stats !== fullStats;
                        if (this.props.trace.name !== '__count') {
                            if (!this.props.trace.continuous) {
                                tooltip += ', mode: ' + stats.value + (showFull ? ' (' + fullStats.value + ')' : '');
                            } else {
                                tooltip += ', z-score: ' + stripTrailingZeros(numberFormat2f(stats.value)) + (showFull ? ' (' + stripTrailingZeros(numberFormat2f(fullStats.value)) + ')' : '');
                            }
                        }
                        tooltip += ', # spots: ' + intFormat(stats.n) + (showFull ? ' / ' + intFormat(fullStats.n) : '');
                    }
                }
                this.props.setTooltip(tooltip);
            } else {
                this.props.setTooltip('');
            }
        });
        containerElement.addEventListener('mouseleave', (e) => {
            this.props.setTooltip('');
        });
        containerElement.addEventListener('click', (e) => {
            if (e.target.nodeName === 'path') {
                const categoryToIndices = this.props.trace.categoryToIndices;
                let category = e.target.id;
                category = category.replaceAll('_', ' '); // FIXME
                const indices = categoryToIndices[category];
                if (indices && indices.length > 0) {
                    this.props.onSelected({
                        name: getEmbeddingKey(this.props.trace.embedding),
                        clear: !e.metaKey && !e.ctrlKey,
                        value: {basis: this.props.trace.embedding, indices: new Set(indices), id: category}
                    });
                } else {
                    this.props.onSelected({name: getEmbeddingKey(this.props.trace.embedding)});
                }
            } else {
                this.props.onSelected({name: getEmbeddingKey(this.props.trace.embedding)});
            }
        });
    }


    componentDidUpdate(prevProps, prevState, snapshot) {
        this.props.setTooltip('');
        if (this.props.trace.source !== prevProps.trace.source) {
            select(prevProps.trace.source).on(".zoom", null);
            this.updateSvg();
        } else if (this.props.chartSize.width !== prevProps.chartSize.width || this.props.chartSize.height !== prevProps.chartSize.height) {
            const svg = this.props.trace.source;
            svg.setAttribute('width', this.props.chartSize.width);
            svg.setAttribute('height', this.props.chartSize.height);
            select(svg).call(this.d3Zoom.transform, zoomIdentity.scale(1));
        }
    }

    render() {
        return <>
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
            </div>

            <div style={{
                display: 'inline-block'
                // width: this.props.chartSize.width,
                // height: this.props.chartSize.height
            }}
                 ref={this.containerElementRef}>
            </div>
        </>;
    }
}


export default withStyles(styles)(MetaEmbedding);

import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import Link from '@material-ui/core/Link';
import Typography from '@material-ui/core/Typography';
import React from 'react';
import ScatterChartThree from './ScatterChartThree';


class GalleryImage extends React.PureComponent {
    constructor(props) {
        super(props);
        this.state = {url: null};
    }


    drawThree() {
        const {traceInfo, scatterGL, containerElement, markerOpacity} = this.props;
        ScatterChartThree.snapshot(scatterGL, traceInfo, markerOpacity);
        const canvas = containerElement.querySelector('canvas');
        // const _this = this;
        // const copy = document.createElement('canvas');
        // let size = {width: 400, height: 400};
        // copy.width = size.width * window.devicePixelRatio;
        // copy.style.width = size.width + 'px';
        // copy.height = size.height * window.devicePixelRatio;
        // copy.style.height = size.height + 'px';
        // const context = copy.getContext('2d');
        // context.scale(window.devicePixelRatio, window.devicePixelRatio);
        // context.drawImage(canvas, 0, 0, size.width, size.height);
        this.setState({url: canvas.toDataURL()});
        // canvas.toBlob(function (blob) {
        //     // let newImg = document.createElement('img');
        //     let url = URL.createObjectURL(blob);
        //     _this.setState({url: url});
        //     // newImg.onload = function () {
        //     //     // no longer need to read the blob so it's revoked
        //     //     URL.revokeObjectURL(url);
        //     // };
        //     //
        //     // newImg.src = url;
        //     // document.body.appendChild(newImg);
        // });

    }


    componentDidMount() {
        this.drawThree();
    }


    componentDidUpdate(prevProps, prevState) {
        this.drawThree();
    }


    onSelect = (event) => {
        event.preventDefault();
        this.props.onSelect(this.props.traceInfo);
    };

    render() {
        let name = this.props.traceInfo.name;
        if (name === '__count') {
            name = 'count';
        }
        return (<Card variant="outlined" style={{display: 'inline-block'}}>
            <CardContent>

                <Typography style={{textAlign: 'center'}} variant="caption" component="h4">
                    <Link href="#"
                          onClick={this.onSelect}>{name}</Link>
                </Typography>
                <img src={this.state.url} style={{width: 400, height: 400}}/>
            </CardContent>
        </Card>);

    }
}

export default GalleryImage;


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
        const url = canvas.toDataURL();
        this.setState({url: url});
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
                <div style={{display: 'inline-block'}}>
                    <img src={this.state.url}/>
                </div>
            </CardContent>

        </Card>);

    }
}

export default GalleryImage;


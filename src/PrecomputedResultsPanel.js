import {Table, TableBody, TableCell, TableHead, TableRow} from '@material-ui/core';
import React from 'react';
import {connect} from 'react-redux';
import {setDistributionPlotOptions} from './actions';


class PrecomputedResultsPanel extends React.PureComponent {

    constructor(props) {
        super(props);
    }

    handleClick = (event, row) => {
        event.stopPropagation();
        // this.props.xx();
    };

    render() {
        const selectedId = null;
        const {dataset} = this.props;
        if (dataset == null || dataset.results == null) {
            return null;
        }
        const columns = dataset.results.columns;
        const rows = dataset.results.rows;
        return <React.Fragment>
            <div style={{maxHeight: 600, overflow: 'auto'}}>
                <Table stickyHeader size={"small"}>
                    <TableHead>
                        <TableRow>
                            {columns.map(column => {
                                return <TableCell key={column.field}>{column.field}</TableCell>;
                            })}
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {rows.map(row => {
                            return <TableRow
                                hover
                                onClick={(event) => this.handleClick(event, row)}
                                role="checkbox"
                                aria-checked={false}
                                tabIndex={-1}
                                key={row.id}
                                selected={row.id === selectedId}
                            >
                                {columns.map(column => {
                                    return <TableCell key={column.field}>{row[column.field]}</TableCell>;
                                })}
                            </TableRow>;
                        })}

                    </TableBody>
                </Table>
            </div>
            <div></div>
        </React.Fragment>;
    }
}


const mapStateToProps = state => {
        return {
            dataset: state.dataset
        };
    }
;
const mapDispatchToProps = (dispatch, ownProps) => {
        return {
            onDistributionPlotOptions: (payload) => {
                dispatch(setDistributionPlotOptions(payload));
            },
        };
    }
;

export default connect(
    mapStateToProps, mapDispatchToProps,
)(PrecomputedResultsPanel);


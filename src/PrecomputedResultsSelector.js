import FormControl from '@material-ui/core/FormControl';
import withStyles from '@material-ui/core/styles/withStyles';
import memoize from 'memoize-one';
import natsort from 'natsort';
import React from 'react';
import {connect} from 'react-redux';
import {setDistributionPlotOptions} from './actions';
import AutocompleteVirtualized from './AutocompleteVirtualized';

const sorter = natsort({insensitive: true});
const getOptions = memoize(
    (rows, field) => {
        const uniqueOptions = new Set();
        rows.forEach(item => {
            uniqueOptions.add(item[field]);
        });
        const options = Array.from(uniqueOptions);
        options.sort(sorter);
        return options;
    }
);
const styles = theme => ({
    formControl: {
        display: 'block',
        minWidth: 200,
        margin: theme.spacing(0, 1)
    }
});
// cirro computed result: {type:['de', 'corr'], params:{}, results:[]}
// for de results is an array of features cluster1_q. cluster_1_fc, cluster_1_percent_expressed, etc.
// facet to pick result, for user submitted jobs: date, name, type
// another facet to filter result (e.g. LFC)
// for custom results
class PrecomputedResultsSelector extends React.PureComponent {

    onChange = (event, value, field) => {
        this.props.dataset.results.facets[field] = value;
        const columns = this.props.dataset.results.columns;
    };

    render() {
        const {dataset, classes} = this.props;
        const results = dataset ? dataset.results : null;
        if (results == null) {
            return null;
        }
        const {columns} = results;

        if (results.facets == null) {
            results.facets = {};
        }
        if (results.allRows == null) {
            results.allRows = results.rows;
        }

        return columns.map(column => {
            return <FormControl key={column.field} className={classes.formControl}>
                <AutocompleteVirtualized
                    label={column.field}
                    multiple={true}
                    options={getOptions(results.allRows, column.field)}
                    value={results.facets[column.field]}
                    onChange={(event, value) => this.onChange(event, value, column.field)}
                />
            </FormControl>;
        });

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


export default withStyles(styles)(connect(
    mapStateToProps, mapDispatchToProps,
)(PrecomputedResultsSelector));

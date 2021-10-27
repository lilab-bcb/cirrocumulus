import {connect} from 'react-redux';
import {find} from 'lodash';
import DotPlotJobResultOptions from './DotPlotJobResultOptions';

function JobResultOptions(props) {
    const {jobResultId, jobResults} = props;
    const jobResult = jobResultId != null ? find(jobResults, item => item.id === jobResultId) : null;
    const jobType = jobResult != null ? jobResult.type : null;
    return <>
        {jobType == 'de' && <DotPlotJobResultOptions jobResult={jobResult}/>}
    </>;
}

const mapStateToProps = state => {
        return {
            jobResultId: state.jobResultId,
            jobResults: state.jobResults
        };
    }
;
const mapDispatchToProps = (dispatch) => {
        return {};
    }
;


export default (connect(
    mapStateToProps, mapDispatchToProps
)(JobResultOptions));
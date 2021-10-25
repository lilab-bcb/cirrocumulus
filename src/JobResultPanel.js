import {connect} from 'react-redux';
import DotPlotJobResultsPanel from './DotPlotJobResultsPanel';
import {find} from 'lodash';
import JobResultsSelector from './JobResultsSelector';

function JobResultPanel(props) {
    const {setTooltip, jobResultId, jobResults} = props;
    const jobResult = jobResultId != null ? find(jobResults, item => item.id === jobResultId) : null;
    const jobType = jobResult != null ? jobResult.type : null;
    return <>
        <JobResultsSelector/>
        {jobType == 'de' && <DotPlotJobResultsPanel setTooltip={setTooltip} jobResult={jobResult}/>}
    </>;
}

const mapStateToProps = state => {
        return {
            jobResultId: state.jobResultId,
            jobResults: state.jobResults,
            tab: state.tab
        };
    }
;
const mapDispatchToProps = (dispatch) => {
        return {};
    }
;


export default (connect(
    mapStateToProps, mapDispatchToProps
)(JobResultPanel));
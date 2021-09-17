import {StyledEngineProvider, ThemeProvider} from '@mui/material/styles';
import {createTheme} from '@mui/material';
import {connect} from 'react-redux';
import App from './App';

const darkTheme = createTheme({palette: {mode: "dark"}});
const lightTheme = createTheme({palette: {mode: "light"}});

function AppWrapper(props) {
    const chartOptions = props.chartOptions;
    const theme = chartOptions.darkMode ? darkTheme : lightTheme;

    return <StyledEngineProvider injectFirst>
        <ThemeProvider theme={theme}>
            <App/>
        </ThemeProvider>
    </StyledEngineProvider>;
}


const mapStateToProps = state => {
    return {
        chartOptions: state.chartOptions
    };
};


export default connect(mapStateToProps, {})(AppWrapper);

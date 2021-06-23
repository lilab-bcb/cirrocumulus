module.exports = {
    preset: "jest-puppeteer",
    testMatch: ["**/tests/*.js"],
    transform: {
        "\\.js$": "react-scripts/config/jest/babelTransform"
    },
    verbose: true,
    testTimeout: 20000
};
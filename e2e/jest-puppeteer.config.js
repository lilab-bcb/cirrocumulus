module.exports = {
    server: {
        command: "python ../cirrocumulus/launch.py ../test-data/pbmc3k_no_raw.h5ad --no-open",
        port: 5000,
        launchTimeout: 10000
    }
}
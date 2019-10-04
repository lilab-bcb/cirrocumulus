class CustomError extends Error {
    constructor(...params) {
        // Pass remaining arguments (including vendor specific ones) to parent constructor
        super(...params);
    }
}

export default CustomError;

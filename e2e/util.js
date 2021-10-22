const gm = require('gm');

module.exports.diffImages = (image1, image2, tolerance) => {
    return new Promise((resolve, reject) => {
        gm.compare(image1, image2, tolerance, function (err, isEqual, equality, raw, path1, path2) {
            resolve(isEqual);
        });
    });
};
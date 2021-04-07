// R chisq.test

import {chiSquare} from './ChiSquare';


it('Agresti(2007) p.39', () => {
    /*
    M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
    dimnames(M) <- list(gender = c("F", "M"),
                        party = c("Democrat","Independent", "Republican"))
    (Xsq <- chisq.test(M, rescale.p=T))
     */
    const counts = [[762, 327, 468], [484, 239, 477]];
    const result = chiSquare(counts);
    expect(result.statistic).toBeCloseTo(30.07015);
    expect(result.p).toBeCloseTo(2.953589e-07);
});
it('test2', () => {
    // counts <- matrix(c(100, 100, 200, 123, 243, 212, 290, 212, 120,234,234,234), nc = 3);
    const result = chiSquare([[100, 243, 120],
        [100, 212, 234],
        [200, 290, 234],
        [123, 212, 234]]);
    expect(result.statistic).toBeCloseTo(57.33124);
    expect(result.p).toBeCloseTo(1.565488e-10);
});
it('test3', () => {
    // counts <- matrix(c(1000, 1000, 1000, 1000, 1000, 1000, 1100, 1000, 1200,1234,1234,1234), nc = 3);
    const result = chiSquare([[1000, 1000, 1200],
        [1000, 1000, 1234],
        [1000, 1100, 1234],
        [1000, 1000, 1234]]);
    expect(result.statistic).toBeCloseTo(4.853864);
    expect(result.p).toBeCloseTo(0.5626884);
});


it('commons-math1', () => {
    // counts <- matrix(c(40, 22, 43, 91, 21, 28, 60, 10, 22), nc = 3);
    const result = chiSquare([[40, 91, 60], [22, 21, 10], [43, 28, 22]]);
    expect(result.statistic).toBeCloseTo(22.709);
    expect(result.p).toBeCloseTo(0.0001447515);

});
it('commons-math2', () => {
    // counts <- matrix(c(10, 15, 30, 40, 60, 90), nc = 3);
    const result = chiSquare([[10, 30, 60], [15, 40, 90]]);
    expect(result.statistic).toBeCloseTo(0.16897);
    expect(result.p).toBeCloseTo(0.919);
});

it('commons-math3', () => {
    // counts <- matrix(c(40, 0, 4, 91, 1, 2, 60, 2, 0), nc = 3);
    const result = chiSquare([[40, 91, 60], [0, 1, 2], [4, 2, 0]]);
    expect(result.statistic).toBeCloseTo(9.6744);
    expect(result.p).toBeCloseTo(0.04628);
});




